// SPDX-License-Identifier: BSD-2-Clause
#include "lut.hpp"
#include "schemelet/rlwe-mp.h"
#include "trees.hpp"
#include <stdexcept>
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
    auto joined = Recompose(s, real, imag);
    const uint32_t decodingLevel = s.depth - s.parameters.levelBudget[1];
    if (joined->GetLevel() > decodingLevel)
        throw std::invalid_argument("Insufficient levels for homomorphic decoding");
    auto decoded = s.cc->EvalHomDecoding(joined, 1, decodingLevel - joined->GetLevel());
    return lbcrypto::SchemeletRLWEMP::ConvertCKKSToRLWE(decoded, s.q);
}
} // namespace large_lut
