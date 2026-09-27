// SPDX-License-Identifier: BSD-2-Clause
#include "digit-decompose.hpp"
#include "schemelet/rlwe-mp.h"
#include "trees.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>

using namespace large_lut;
using namespace lbcrypto;

namespace {
void Require(bool condition, const std::string& message) {
    if (!condition)
        throw std::runtime_error(message);
}

template <class Function> void RequireInvalid(Function&& function) {
    try {
        function();
    } catch (const std::invalid_argument&) {
        return;
    }
    throw std::runtime_error("Expected invalid_argument for malformed ciphertext");
}

std::vector<int64_t> MakeInputs(const Setup& s, uint64_t seed) {
    std::mt19937_64 random(seed);
    std::vector<int64_t> values(s.parameters.ringDim);
    // Exhaust small domains, including every carry boundary.
    if (s.lutSize <= values.size()) {
        for (size_t i = 0; i < values.size(); ++i)
            values[i] = i % s.lutSize;
        std::shuffle(values.begin(), values.end(), random);
        return values;
    } else {
        for (size_t i = 0; i < values.size(); ++i) {
            uint64_t value = 0;
            uint32_t shift = 0;
            for (size_t digit = 0; digit < s.digitWidths.size(); ++digit) {
                uint64_t mask = (1u << s.digitWidths[digit]) - 1;
                value |= ((i * (2 * digit + 1) + seed + digit) & mask) << shift;
                shift += s.digitWidths[digit];
            }
            values[i] = value;
        }
        // Retain random combinations as well as carry boundaries.
        for (size_t i = 0; i < values.size(); i += 8)
            values[i] = random() % s.lutSize;
    }
    size_t index = 0;
    for (uint32_t bit = 0; bit < s.parameters.inputBits; ++bit) {
        const int64_t boundary = int64_t{1} << bit;
        values[index++] = boundary - 1;
        values[index++] = boundary;
        values[index++] = (boundary + 1) % s.lutSize;
    }
    values[index] = s.lutSize - 1;
    return values;
}

std::vector<int64_t> Decode(const Setup& s, const large_lut::Ciphertext& real,
                            const large_lut::Ciphertext& imag) {
    const auto realElements = real->GetElements();
    const auto imagElements = imag->GetElements();
    auto joined = Recompose(s, real, imag);
    Require(real->GetElements() == realElements && imag->GetElements() == imagElements,
            "Recomposition changed its inputs");
    const auto targetLevel = s.depth - s.parameters.levelBudget[1];
    auto decoded = s.cc->EvalHomDecoding(joined, 1, targetLevel - joined->GetLevel());
    return DecryptOutput(s, SchemeletRLWEMP::ConvertCKKSToRLWE(decoded, s.q));
}

size_t CoefficientIndex(size_t slot, uint32_t slots) {
    size_t reversed = 0;
    for (uint32_t width = slots; width > 1; width >>= 1) {
        reversed = (reversed << 1) | (slot & 1);
        slot >>= 1;
    }
    return reversed;
}

double CheckSelectors(const Setup& s, const SplitSelectors& selectors, const std::vector<int64_t>& input) {
    Require(selectors.real.size() == s.digitWidths.size(), "Wrong real digit count");
    Require(selectors.imag.size() == s.digitWidths.size(), "Wrong imaginary digit count");
    double maxError = 0;
    uint32_t shift = 0;
    for (size_t digit = 0; digit < s.digitWidths.size(); ++digit) {
        const uint32_t mask = (1u << s.digitWidths[digit]) - 1;
        Require(selectors.real[digit].size() == mask, "Wrong real selector count");
        Require(selectors.imag[digit].size() == mask, "Wrong imaginary selector count");
        for (uint32_t selector = 0; selector < mask; ++selector) {
            for (size_t half = 0; half < 2; ++half) {
                const auto& ct =
                    half == 0 ? selectors.real[digit][selector] : selectors.imag[digit][selector];
                Plaintext plaintext;
                s.cc->Decrypt(s.secretKey, ct, &plaintext);
                plaintext->SetLength(s.slots);
                auto values = plaintext->GetCKKSPackedValue();
                for (size_t slot = 0; slot < s.slots; ++slot) {
                    const auto expected =
                        double(((uint64_t(input[half * s.slots + CoefficientIndex(slot, s.slots)]) >> shift) &
                                mask) == selector);
                    const double error = std::abs(values[slot] - expected);
                    Require(std::isfinite(error) && error < 1e-4,
                            "Selector error at digit " + std::to_string(digit) + ", selector " +
                                std::to_string(selector) + ", half " + std::to_string(half) + ", slot " +
                                std::to_string(slot) + ": " + std::to_string(error));
                    maxError = std::max(maxError, error);
                }
            }
        }

        // Reconstruct the omitted selector to check partition of unity. Cleaning
        // is unnecessary for this unweighted decoding check.
        auto realSum = selectors.real[digit][0]->Clone();
        auto imagSum = selectors.imag[digit][0]->Clone();
        for (uint32_t selector = 1; selector < mask; ++selector) {
            s.cc->EvalAddInPlace(realSum, selectors.real[digit][selector]);
            s.cc->EvalAddInPlace(imagSum, selectors.imag[digit][selector]);
        }
        auto lastReal = s.cc->EvalSub(1.0, realSum);
        auto lastImag = s.cc->EvalSub(1.0, imagSum);
        auto decoded = Decode(s, lastReal, lastImag);
        for (size_t i = 0; i < input.size(); ++i)
            Require(decoded[i] == int64_t(((uint64_t(input[i]) >> shift) & mask) == mask),
                    "Omitted selector decoding mismatch");

        // Recompose and decode a digit from its selectors, checking scale and
        // packing independently of direct CKKS slot decryption.
        auto realDigit = s.cc->GetScheme()->MultByInteger(lastReal, mask);
        auto imagDigit = s.cc->GetScheme()->MultByInteger(lastImag, mask);
        for (uint32_t selector = 1; selector < mask; ++selector) {
            s.cc->EvalAddInPlace(realDigit,
                                 s.cc->GetScheme()->MultByInteger(selectors.real[digit][selector], selector));
            s.cc->EvalAddInPlace(imagDigit,
                                 s.cc->GetScheme()->MultByInteger(selectors.imag[digit][selector], selector));
        }
        decoded = Decode(s, realDigit, imagDigit);
        for (size_t i = 0; i < input.size(); ++i)
            Require(decoded[i] == int64_t((uint64_t(input[i]) >> shift) & mask), "Digit decoding mismatch");
        shift += s.digitWidths[digit];
    }
    return maxError;
}

void TestCase(uint32_t inputBits, uint32_t digitBits, uint32_t h = 192, bool repeat = false,
              bool earlyDrops = false, LutKind kind = LutKind::Cleartext) {
    std::cout << "Input " << inputBits << ", digit " << digitBits << ", h " << h << ", early drops "
              << earlyDrops << std::flush;
    Parameters parameters;
    parameters.inputBits = inputBits;
    parameters.digitBits = digitBits;
    parameters.h = h;
    parameters.lutKind = kind;
    parameters.toy = true;
    parameters.ringDim = 256;
    auto s = CreateSetup(parameters);
    auto evaluator = s;
    evaluator.secretKey.reset();
    for (uint64_t seed = 1; seed <= (repeat ? 2u : 1u); ++seed) {
        auto values = MakeInputs(s, seed);
        auto input = EncryptInput(s, values);
        auto unchanged = input;
        const auto drops = earlyDrops ? large_lut::detail::SelectorLevelDrops(s) : std::vector<uint32_t>{};
        auto selectors = Decompose(evaluator, input, drops);
        const auto baseLevel = s.depth - s.computationDepth - s.parameters.levelBudget[1];
        for (const auto* half : {&selectors.real, &selectors.imag})
            for (size_t digit = 0; digit < half->size(); ++digit)
                for (const auto& ct : (*half)[digit])
                    Require(ct->GetLevel() == baseLevel + (earlyDrops ? drops[digit] : 0) &&
                                ct->GetNoiseScaleDeg() == 1,
                            "Unexpected selector level or scale degree");
        Require(input == unchanged, "Decomposition changed its input");
        const double error = CheckSelectors(s, selectors, values);
        std::cout << ", error " << error << std::flush;
        if (repeat) {
            RequireInvalid([&] { Decompose(evaluator, {}); });
            const std::vector<uint32_t> wrongSize(s.digitWidths.size() + 1);
            RequireInvalid([&] { Decompose(evaluator, input, wrongSize); });
            const std::vector<uint32_t> tooDeep(s.digitWidths.size(), s.computationDepth + 1);
            RequireInvalid([&] { Decompose(evaluator, input, tooDeep); });
            auto malformed = input;
            malformed[0].SwitchModulus(s.q >> 1, 1, 0, 0);
            RequireInvalid([&] { Decompose(evaluator, malformed); });
        }
    }
    std::cout << '\n';
    s.cc->ClearEvalMultKeys();
    s.cc->ClearEvalAutomorphismKeys();
    CryptoContextFactory<DCRTPoly>::ReleaseAllContexts();
}
} // namespace

int main() {
    try {
        TestCase(1, 8);
        for (uint32_t bits = 1; bits <= 8; ++bits) {
            TestCase(2 * bits, bits, 192, bits <= 2);
            if (bits > 1)
                TestCase(2 * bits + 1, bits);
        }
        for (uint32_t top = 2; top < 8; ++top)
            TestCase(16 + top, 8);
        for (uint32_t bits = 1; bits <= 8; ++bits)
            TestCase(32, bits);
        TestCase(15, 4);
        TestCase(9, 3, 192, true, true);
        TestCase(17, 4, 192, false, true);
        TestCase(9, 8, 192, false, true, LutKind::Ciphertext);
        TestCase(5, 2, 192, false, true, LutKind::Plaintext);
        TestCase(32, 1, 192, false, true);
        TestCase(31, 1, 192, false, true);
        TestCase(13, 4, 32);
        TestCase(13, 4, 256);
    } catch (const std::exception& error) {
        std::cerr << "\n" << error.what() << '\n';
        return 1;
    }
    std::cout << "Digit decomposition tests passed\n";
}
