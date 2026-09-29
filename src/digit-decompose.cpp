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

void DropPowers(const Setup& s, const SplitPowers& powers, uint32_t levels) {
    if (levels == 0)
        return;
    for (const auto& half : powers) {
        for (auto* row : {&half->powersRe, &half->powers2Re})
            for (auto& power : *row)
                if (power)
                    power = s.cc->LevelReduce(power, nullptr, levels);
        if (half->power2km1Re)
            half->power2km1Re = s.cc->LevelReduce(half->power2km1Re, nullptr, levels);
    }
}

Ciphertext EvaluateHalf(const Setup& s, const Powers& powers, const Coefficients& coefficients) {
    if (s.digitBits == 1) {
        // Binary precomputation gives cos^2(pi*x/2). Clone before applying
        // the affine map so the other functions can reuse the same powers.
        auto result = powers->powersRe.front()->Clone();
        return coefficients[1].real() < 0 ? s.cc->EvalSub(1.0, result) : result;
    }
    auto result = s.cc->EvalPolyWithPrecomp(powers, coefficients);
    auto conjugate = FHECKKSRNS::Conjugate(result, s.cc->GetEvalAutomorphismKeyMap(result->GetKeyTag()));
    // The Hermite coefficients already include the factor 1/2.
    s.cc->EvalAddInPlace(result, conjugate);
    return result;
}
} // namespace

Ciphertext Recompose(const Setup& s, const Ciphertext& real, const Ciphertext& imag) {
    auto result = s.cc->GetScheme()->MultByMonomial(imag, s.cc->GetCyclotomicOrder() / 4);
    s.cc->EvalAddInPlace(result, real);
    return result;
}

SplitSelectors Decompose(const Setup& s, const RLWECiphertext& input, std::span<const uint32_t> levelDrops) {
    if (!levelDrops.empty() && levelDrops.size() != s.digitWidths.size())
        throw std::invalid_argument("Expected one level drop per digit");
    for (auto drop : levelDrops)
        if (drop > s.computationDepth)
            throw std::invalid_argument("Digit level drop exceeds the computation depth");
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

        const uint32_t drop = levelDrops.empty() ? 0 : levelDrops[digit];
        DropPowers(s, powers, drop);
        const bool hasNext = digit + 1 < s.digitWidths.size();
        SplitCiphertext extracted;

        const uint32_t count = (1u << s.digitWidths[digit]) - 1;
        const uint32_t stride = 1u << (s.digitBits - s.digitWidths[digit]);
        result.real[digit].resize(count);
        result.imag[digit].resize(count);
        const uint32_t functions = count + hasNext;
#pragma omp parallel for schedule(static)
        for (uint32_t task = 0; task < 2 * functions; ++task) {
            const uint32_t function = task / 2;
            const uint32_t half = task % 2;
            const auto& coefficients =
                function == count ? s.modCoefficients : s.selectorCoefficients[function * stride];
            auto value = EvaluateHalf(s, powers[half], coefficients);
            if (function == count)
                extracted[half] = std::move(value);
            else
                (half == 0 ? result.real : result.imag)[digit][function] = std::move(value);
        }

        if (hasNext) {
            auto joined = Recompose(s, extracted[0], extracted[1]);
            // Decoding was set up for the original input modulus. Restore the
            // current residual's scale, including the digit polynomial's scale.
            auto decoded =
                s.cc->EvalHomDecoding(joined, uint64_t(s.modScale) << removedBits, s.computationDepth - drop);
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
