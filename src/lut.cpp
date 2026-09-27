// SPDX-License-Identifier: BSD-2-Clause
#include "lut.hpp"
#include "schemelet/rlwe-mp.h"
#include "trees.hpp"
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace large_lut {
namespace {
Ciphertext Clean(const Setup& s, const Ciphertext& input) {
    auto square = s.cc->EvalSquare(input);
    s.cc->ModReduceInPlace(square);
    auto affine = s.cc->EvalSub(3.0, s.cc->GetScheme()->MultByInteger(input, 2));
    auto result = s.cc->EvalMult(square, affine);
    s.cc->ModReduceInPlace(result);
    return result;
}
RLWECiphertext DecodeResult(const Setup& s, const Ciphertext& real, const Ciphertext& imag) {
    auto joined = Recompose(s, real, imag);
    const uint32_t decodingLevel = s.depth - s.parameters.levelBudget[1];
    if (joined->GetLevel() > decodingLevel)
        throw std::invalid_argument("Insufficient levels for homomorphic decoding");
    auto decoded = s.cc->EvalHomDecoding(joined, 1, decodingLevel - joined->GetLevel());
    return lbcrypto::SchemeletRLWEMP::ConvertCKKSToRLWE(decoded, s.q);
}

void CheckPlaintext(const Setup& s, const lbcrypto::Plaintext& entry) {
    if (!entry || entry->GetEncodingType() != lbcrypto::CKKS_PACKED_ENCODING ||
        entry->GetLevel() != detail::LutEntryLevel(s) || entry->GetNoiseScaleDeg() != 1 ||
        entry->GetSlots() != s.slots)
        throw std::invalid_argument("Invalid plaintext LUT entry metadata");
    const auto& poly = entry->GetElement<lbcrypto::DCRTPoly>();
    const auto& params = s.cc->GetCryptoParameters()->GetElementParams()->GetParams();
    if (poly.GetRingDimension() != s.parameters.ringDim ||
        poly.GetNumOfElements() != params.size() - entry->GetLevel())
        throw std::invalid_argument("Invalid plaintext LUT entry parameters");
    for (size_t i = 0; i < poly.GetNumOfElements(); ++i)
        if (*poly.GetParams()->GetParams()[i] != *params[i])
            throw std::invalid_argument("Plaintext LUT entry has a different modulus chain");
    auto crypto = std::dynamic_pointer_cast<lbcrypto::CryptoParametersCKKSRNS>(s.cc->GetCryptoParameters());
    if (entry->GetScalingFactor() != crypto->GetScalingFactorReal(entry->GetLevel()))
        throw std::invalid_argument("Invalid plaintext LUT entry scale");
}

template <class Lut> void CheckEncodedLut(const Setup& s, const Lut& lut) {
    if (lut.real.size() != s.lutSize || lut.imag.size() != s.lutSize)
        throw std::invalid_argument("Both LUT halves must contain 2^inputBits entries");
    for (const auto* half : {&lut.real, &lut.imag}) {
        for (const auto& entry : *half) {
            if constexpr (std::is_same_v<Lut, PlaintextLut>) {
                CheckPlaintext(s, entry);
            } else {
                if (!entry || entry->GetCryptoContext() != s.cc ||
                    entry->GetKeyTag() != s.publicKey->GetKeyTag() ||
                    entry->GetLevel() != detail::LutEntryLevel(s) || entry->GetNoiseScaleDeg() != 1 ||
                    entry->GetSlots() != s.slots || entry->GetElements().size() != 2)
                    throw std::invalid_argument("Invalid ciphertext LUT entry");
                auto crypto =
                    std::dynamic_pointer_cast<lbcrypto::CryptoParametersCKKSRNS>(s.cc->GetCryptoParameters());
                if (entry->GetScalingFactor() != crypto->GetScalingFactorReal(entry->GetLevel()))
                    throw std::invalid_argument("Invalid ciphertext LUT entry scale");
            }
        }
    }
}

} // namespace

SplitSelectors PrepareSelectors(const Setup& s, SplitSelectors selectors) {
    for (const auto* matrix : {&selectors.real, &selectors.imag}) {
        if (matrix->size() != s.digitWidths.size())
            throw std::invalid_argument("Unexpected selector digit count");
        for (size_t digit = 0; digit < matrix->size(); ++digit) {
            const auto& row = (*matrix)[digit];
            if (row.size() != (1u << s.digitWidths[digit]) - 1)
                throw std::invalid_argument("Expected all but one selector per digit");
            for (const auto& ct : row)
                if (!ct || ct->GetCryptoContext() != s.cc || ct->GetKeyTag() != s.publicKey->GetKeyTag())
                    throw std::invalid_argument("Invalid selector ciphertext");
        }
    }
    for (auto* matrix : {&selectors.real, &selectors.imag}) {
        for (auto& row : *matrix) {
#pragma omp parallel for
            for (size_t i = 0; i < row.size(); ++i)
                for (uint32_t iteration = 0; iteration < s.cleaning; ++iteration)
                    row[i] = Clean(s, row[i]);
            auto sum = row[0]->Clone();
            for (size_t i = 1; i < row.size(); ++i)
                s.cc->EvalAddInPlace(sum, row[i]);
            row.push_back(s.cc->EvalSub(1.0, sum));
        }
    }
    return selectors;
}

RLWECiphertext Evaluate(const Setup& s, const RLWECiphertext& input, std::span<const int64_t> lut) {
    if (s.parameters.lutKind != LutKind::Cleartext)
        throw std::invalid_argument("Integer LUT entries require a cleartext setup");
    if (lut.size() != s.lutSize)
        throw std::invalid_argument("LUT size must be 2^inputBits");
    for (auto value : lut)
        if (value < 0 || uint64_t(value) >= s.lutSize)
            throw std::invalid_argument("LUT output outside [0, 2^inputBits)");

    auto selectors = PrepareSelectors(s, Decompose(s, input));
    auto real = EvaluateCleartextTree(s, std::move(selectors.real), lut);
    auto imag = EvaluateCleartextTree(s, std::move(selectors.imag), lut);
    return DecodeResult(s, real, imag);
}
PlaintextLut EncodeLut(const Setup& s, const LutValues& values) {
    if (s.parameters.lutKind == LutKind::Cleartext || !values)
        throw std::invalid_argument("Encoding requires a plaintext or ciphertext setup and a value callback");
    PlaintextLut result;
    result.real.reserve(s.lutSize);
    result.imag.reserve(s.lutSize);
    std::vector<uint32_t> order(s.slots);
    for (uint32_t slot = 0; slot < s.slots; ++slot) {
        auto index = slot;
        for (uint32_t width = s.slots; width > 1; width >>= 1) {
            order[slot] = (order[slot] << 1) | (index & 1);
            index >>= 1;
        }
    }
    for (uint64_t index = 0; index < s.lutSize; ++index) {
        for (uint32_t half = 0; half < 2; ++half) {
            std::vector<double> packed(s.slots);
            for (uint32_t slot = 0; slot < s.slots; ++slot) {
                const int64_t value = values(index, half * s.slots + order[slot]);
                if (value < 0 || uint64_t(value) >= s.lutSize)
                    throw std::invalid_argument("LUT output outside [0, 2^inputBits)");
                packed[slot] = value;
            }
            auto entry = s.cc->MakeCKKSPackedPlaintext(packed, 1, detail::LutEntryLevel(s));
            (half == 0 ? result.real : result.imag).push_back(std::move(entry));
        }
    }
    return result;
}

PlaintextLut EncodeLut(const Setup& s, std::span<const int64_t> values) {
    if (values.size() != s.lutSize)
        throw std::invalid_argument("LUT size must be 2^inputBits");
    return EncodeLut(s, [&](uint64_t index, uint32_t) { return values[index]; });
}

CiphertextLut EncryptLut(const Setup& s, const PlaintextLut& lut) {
    if (s.parameters.lutKind != LutKind::Ciphertext)
        throw std::invalid_argument("Encryption requires a ciphertext LUT setup");
    CheckEncodedLut(s, lut);
    CiphertextLut result;
    for (const auto* half : {&lut.real, &lut.imag}) {
        auto& output = half == &lut.real ? result.real : result.imag;
        output.reserve(s.lutSize);
        for (const auto& entry : *half)
            output.push_back(s.cc->Encrypt(s.publicKey, entry));
    }
    return result;
}

RLWECiphertext Evaluate(const Setup& s, const RLWECiphertext& input, const PlaintextLut& lut) {
    if (s.parameters.lutKind != LutKind::Plaintext)
        throw std::invalid_argument("Plaintext entries require a plaintext LUT setup");
    CheckEncodedLut(s, lut);
    auto selectors = PrepareSelectors(s, Decompose(s, input));
    auto real = detail::EvaluatePlaintextTree(s, std::move(selectors.real), lut.real);
    auto imag = detail::EvaluatePlaintextTree(s, std::move(selectors.imag), lut.imag);
    return DecodeResult(s, real, imag);
}

RLWECiphertext Evaluate(const Setup& s, const RLWECiphertext& input, const CiphertextLut& lut) {
    if (s.parameters.lutKind != LutKind::Ciphertext)
        throw std::invalid_argument("Ciphertext entries require a ciphertext LUT setup");
    CheckEncodedLut(s, lut);
    auto selectors = PrepareSelectors(s, Decompose(s, input));
    auto real = detail::EvaluateCiphertextTree(s, std::move(selectors.real), lut.real);
    auto imag = detail::EvaluateCiphertextTree(s, std::move(selectors.imag), lut.imag);
    return DecodeResult(s, real, imag);
}
} // namespace large_lut
