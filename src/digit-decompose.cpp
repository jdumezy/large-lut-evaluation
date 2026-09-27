// SPDX-License-Identifier: BSD-2-Clause
#include "digit-decompose.hpp"
#include "scheme/ckksrns/ckksrns-fhe.h"
#include "schemelet/rlwe-mp.h"
#include <array>
#include <stdexcept>
#include <utility>

namespace large_lut {
namespace {
using namespace lbcrypto;
using Powers = std::shared_ptr<seriesPowers<DCRTPoly>>;
using SplitPowers = std::array<Powers, 2>;
using SplitCiphertext = std::array<Ciphertext, 2>;

SplitPowers SeparatePowers(const Powers& powers) {
    if (powers->powers2Re.empty())
        return {std::make_shared<seriesPowers<DCRTPoly>>(powers->powersRe),
                std::make_shared<seriesPowers<DCRTPoly>>(powers->powersIm)};
    return {std::make_shared<seriesPowers<DCRTPoly>>(powers->powersRe, powers->powers2Re, powers->power2km1Re,
                                                     powers->k, powers->m),
            std::make_shared<seriesPowers<DCRTPoly>>(powers->powersIm, powers->powers2Im, powers->power2km1Im,
                                                     powers->k, powers->m)};
}

SplitCiphertext EvaluateSplit(const Setup& s, const SplitPowers& powers, const Coefficients& coefficients) {
    SplitCiphertext result;
    for (size_t half = 0; half < result.size(); ++half) {
        if (s.digitBits == 1) {
            // Binary precomputation gives cos^2(pi*x/2). Clone before applying
            // the affine map so the other functions can reuse the same powers.
            result[half] = powers[half]->powersRe.front()->Clone();
            if (coefficients[1].real() < 0)
                result[half] = s.cc->EvalSub(1.0, result[half]);
        } else {
            auto localPowers = powers[half];
            if (coefficients.size() < 6) {
                // OpenFHE's linear evaluator scales its input powers in place.
                auto copies = powers[half]->powersRe;
                for (auto& power : copies)
                    power = power->Clone();
                localPowers = std::make_shared<seriesPowers<DCRTPoly>>(copies);
            }
            result[half] = s.cc->EvalPolyWithPrecomp(localPowers, coefficients);
            auto conjugate = FHECKKSRNS::Conjugate(
                result[half], s.cc->GetEvalAutomorphismKeyMap(result[half]->GetKeyTag()));
            // The Hermite coefficients already include the factor 1/2.
            s.cc->EvalAddInPlace(result[half], conjugate);
        }
    }
    return result;
}
} // namespace

Ciphertext Recompose(const Setup& s, const Ciphertext& real, const Ciphertext& imag) {
    auto result = s.cc->GetScheme()->MultByMonomial(imag, s.cc->GetCyclotomicOrder() / 4);
    s.cc->EvalAddInPlace(result, real);
    return result;
}

SplitSelectors Decompose(const Setup& s, const RLWECiphertext& input) {
    if (input.size() != 2)
        throw std::invalid_argument("An RLWE ciphertext must have two polynomials");
    for (const auto& poly : input)
        if (poly.GetRingDimension() != s.parameters.ringDim || poly.GetModulus() != s.inputQ ||
            poly.GetFormat() != COEFFICIENT)
            throw std::invalid_argument("Unexpected RLWE ring dimension, modulus or format");

    SplitSelectors result{SelectorMatrix(s.digitWidths.size()), SelectorMatrix(s.digitWidths.size())};
    auto residual = input;
    auto currentQ = s.inputQ;
    uint32_t removedBits = 0;
    for (size_t digit = 0; digit < s.digitWidths.size(); ++digit) {
        auto reduced = residual;
        for (auto& poly : reduced)
            poly.SwitchModulus(s.q, 1, 0, 0);
        auto ckks = SchemeletRLWEMP::ConvertRLWEToCKKS(*s.cc, reduced, s.publicKey, s.q, s.slots, s.depth);
        auto powers = SeparatePowers(s.cc->EvalMVBPrecompute(ckks, s.modCoefficients, s.digitBits, s.qPrime));

        const bool hasNext = digit + 1 < s.digitWidths.size();
        SplitCiphertext extracted;
        if (hasNext)
            extracted = EvaluateSplit(s, powers, s.modCoefficients);

        const uint32_t count = (1u << s.digitWidths[digit]) - 1;
        const uint32_t stride = 1u << (s.digitBits - s.digitWidths[digit]);
        result.real[digit].resize(count);
        result.imag[digit].resize(count);
#pragma omp parallel for
        for (uint32_t selector = 0; selector < count; ++selector) {
            auto split = EvaluateSplit(s, powers, s.selectorCoefficients[selector * stride]);
            result.real[digit][selector] = std::move(split[0]);
            result.imag[digit][selector] = std::move(split[1]);
        }

        if (hasNext) {
            auto joined = Recompose(s, extracted[0], extracted[1]);
            // Decoding was set up for the original input modulus. Restore the
            // current residual's scale, including the digit polynomial's scale.
            auto decoded =
                s.cc->EvalHomDecoding(joined, uint64_t(s.modScale) << removedBits, s.computationDepth);
            auto lowDigit = SchemeletRLWEMP::ConvertCKKSToRLWE(decoded, currentQ);
            for (size_t component = 0; component < residual.size(); ++component)
                residual[component] -= lowDigit[component];

            // Reducing by the next width leaves q for the final extraction.
            // A partial top digit is embedded at multiples of 2^(b-w), matching
            // the stride used to select its polynomials above.
            const uint32_t nextWidth = s.digitWidths[digit + 1];
            auto oldQ = currentQ;
            currentQ >>= nextWidth;
            for (auto& poly : residual) {
                poly = poly.MultiplyAndRound(currentQ, oldQ);
                poly.SwitchModulus(currentQ, 1, 0, 0);
            }
            removedBits += nextWidth;
        }
    }
    return result;
}
} // namespace large_lut
