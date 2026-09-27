// SPDX-License-Identifier: BSD-2-Clause
#include "lut.hpp"
#include "trees.hpp"
#include <algorithm>
#include <bit>
#include <cmath>
#include <iostream>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>

using namespace large_lut;
using lbcrypto::CryptoContextFactory;
using lbcrypto::DCRTPoly;
using lbcrypto::Plaintext;

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
    throw std::runtime_error("Expected invalid_argument");
}

Parameters ToyParameters(uint32_t inputBits, uint32_t digitBits, int cleaning = -1) {
    Parameters p;
    p.toy = true;
    p.ringDim = 256;
    p.inputBits = inputBits;
    p.digitBits = digitBits;
    p.cleaning = cleaning;
    return p;
}

void Release(const Setup& s) {
    s.cc->ClearEvalMultKeys();
    s.cc->ClearEvalAutomorphismKeys();
    CryptoContextFactory<DCRTPoly>::ReleaseAllContexts();
}

std::vector<std::complex<double>> Slots(const Setup& s, const Ciphertext& ct) {
    Plaintext plaintext;
    s.cc->Decrypt(s.secretKey, ct, &plaintext);
    plaintext->SetLength(s.slots);
    return plaintext->GetCKKSPackedValue();
}

void TestSchedules() {
    const std::vector<std::vector<uint32_t>> expected = {
        {0},
        {0, 0},
        {0, 0, 1},
        {0, 0, 0, 0},
        {0, 0, 1, 1, 1},
        {0, 0, 1, 0, 0, 1},
        {0, 0, 0, 0, 0, 0, 1},
        {0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 2, 1, 1, 1, 1},
        {0, 0, 0, 0, 2, 0, 0, 0, 0, 2},
        {0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 2},
        {0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1},
        {0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1},
        {0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 3, 1, 1, 1, 1, 1, 1, 1, 1},
        {0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 3},
        {0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 3},
        {0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2},
        {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2},
        {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2},
        {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2},
        {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1}};
    for (size_t i = 0; i < expected.size(); ++i)
        Require(TreeLevelDrops(i + 1) == expected[i], "Incorrect tree level schedule");
    RequireInvalid([] { TreeLevelDrops(0); });
    RequireInvalid([] { TreeLevelDrops(33); });
}

void TestCleaning() {
    auto s = CreateSetup(ToyParameters(8, 3, 2));
    const uint32_t level = s.depth - s.computationDepth - s.parameters.levelBudget[1];
    SplitSelectors raw{SelectorMatrix(s.digitWidths.size()), SelectorMatrix(s.digitWidths.size())};
    for (size_t half = 0; half < 2; ++half) {
        auto& matrix = half == 0 ? raw.real : raw.imag;
        for (size_t digit = 0; digit < matrix.size(); ++digit) {
            const uint32_t base = 1u << s.digitWidths[digit];
            for (uint32_t selector = 0; selector + 1 < base; ++selector) {
                std::vector<double> values(s.slots);
                for (size_t slot = 0; slot < values.size(); ++slot)
                    values[slot] = ((slot + half + digit) % base == selector) ? 0.98 : 0.01;
                matrix[digit].push_back(
                    s.cc->Encrypt(s.publicKey, s.cc->MakeCKKSPackedPlaintext(values, 1, level)));
            }
        }
    }
    auto original = raw.real[0][0]->GetElements();
    auto clean = PrepareSelectors(s, raw);
    Require(raw.real[0][0]->GetElements() == original, "Cleaning mutated its input");
    for (size_t half = 0; half < 2; ++half) {
        const auto& matrix = half == 0 ? clean.real : clean.imag;
        for (size_t digit = 0; digit < matrix.size(); ++digit) {
            const uint32_t base = 1u << s.digitWidths[digit];
            Require(matrix[digit].size() == base, "Missing complement selector");
            std::vector<double> sum(s.slots);
            for (uint32_t selector = 0; selector < base; ++selector) {
                const auto& ct = matrix[digit][selector];
                Require(ct->GetLevel() == level + 2 * s.cleaning, "Incorrect cleaning depth");
                auto values = Slots(s, ct);
                for (size_t slot = 0; slot < values.size(); ++slot) {
                    double expected = 1 - sum[slot];
                    if (selector + 1 < base) {
                        expected = ((slot + half + digit) % base == selector) ? 0.98 : 0.01;
                        for (uint32_t iteration = 0; iteration < s.cleaning; ++iteration)
                            expected = expected * expected * (3 - 2 * expected);
                        sum[slot] += expected;
                    }
                    Require(std::abs(values[slot] - expected) < 1e-7,
                            "Cleaning or complement order is incorrect");
                }
            }
        }
    }
    RequireInvalid([&] { PrepareSelectors(s, clean); });
    RequireInvalid([&] { PrepareSelectors(s, {}); });
    Release(s);
    std::cout << "Cleaning and complement checks passed\n";
}

void TestTree(uint32_t digits) {
    auto s = CreateSetup(ToyParameters(digits, 1, 0));
    const uint32_t level = s.depth - s.computationDepth - s.parameters.levelBudget[1];
    SelectorMatrix selectors(digits);
    std::vector<size_t> indices(s.slots);
    size_t entries = 1;
    // Narrow rows exercise every multiplication schedule through 32 digits
    // without allocating an exponential LUT. Three rows have two choices.
    for (uint32_t digit = 0; digit < digits; ++digit) {
        const size_t width = (digit == 0 || digit == digits / 2 || digit + 1 == digits) ? 2 : 1;
        for (size_t slot = 0; slot < s.slots; ++slot)
            indices[slot] += ((slot >> (digit % 7)) % width) * entries;
        for (size_t selector = 0; selector < width; ++selector) {
            std::vector<double> values(s.slots);
            for (size_t slot = 0; slot < s.slots; ++slot)
                values[slot] = ((slot >> (digit % 7)) % width) == selector;
            selectors[digit].push_back(
                s.cc->Encrypt(s.publicKey, s.cc->MakeCKKSPackedPlaintext(values, 1, level)));
        }
        entries *= width;
    }
    std::vector<int64_t> lut(entries);
    for (size_t i = 0; i < entries; ++i)
        lut[i] = 17 + 37 * i;
    const auto original = selectors.back()[0]->GetElements();
    auto result = EvaluateCleartextTree(s, selectors, lut);
    Require(selectors.back()[0]->GetElements() == original, "Tree mutated its selectors");
    Require(result->GetLevel() == level + std::bit_width(digits - 1), "Tree depth mismatch");
    Require(result->GetElements().size() == 2 && result->GetNoiseScaleDeg() == 1,
            "Tree did not relinearize and rescale");
    auto values = Slots(s, result);
    for (size_t slot = 0; slot < s.slots; ++slot)
        Require(std::abs(values[slot] - double(lut[indices[slot]])) < 1e-6, "Tree lookup mismatch");
    RequireInvalid([&] { EvaluateCleartextTree(s, selectors, {}); });
    Release(s);
}

void TestLookup(uint32_t inputBits, uint32_t digitBits, int cleaning = -1, bool special = false) {
    std::cout << "LUT " << inputBits << " bits, digit " << digitBits << ", cleaning " << cleaning
              << std::flush;
    auto s = CreateSetup(ToyParameters(inputBits, digitBits, cleaning));
    auto evaluator = s;
    evaluator.secretKey.reset();
    std::mt19937_64 random(101 + inputBits * 10 + digitBits);
    std::vector<int64_t> lut(s.lutSize);
    for (auto& value : lut)
        value = random() % s.lutSize;
    std::vector<int64_t> values(s.parameters.ringDim);
    for (size_t i = 0; i < values.size(); ++i)
        values[i] = s.lutSize <= values.size() ? i % s.lutSize : random() % s.lutSize;
    if (s.lutSize > values.size()) {
        size_t i = 0;
        for (uint32_t bit = 0; bit < inputBits; ++bit) {
            values[i++] = (int64_t{1} << bit) - 1;
            values[i++] = int64_t{1} << bit;
        }
    }
    values.back() = s.lutSize - 1;
    auto input = EncryptInput(s, values);
    const auto unchanged = input;
    auto check = [&] {
        auto output = Evaluate(evaluator, input, lut);
        Require(input == unchanged, "LUT evaluation mutated its input");
        Require(output.size() == 2 && output[0].GetModulus() == s.q && output[1].GetModulus() == s.q,
                "Wrong RLWE output format");
        auto actual = DecryptOutput(s, output);
        for (size_t i = 0; i < values.size(); ++i)
            Require(actual[i] == lut[values[i]], "Lookup mismatch at slot " + std::to_string(i) + ": " +
                                                     std::to_string(actual[i]) +
                                                     " != " + std::to_string(lut[values[i]]));
    };
    check();
    if (special) {
        std::iota(lut.begin(), lut.end(), 0);
        check();
        std::fill(lut.begin(), lut.end(), 0);
        check();
        std::fill(lut.begin(), lut.end(), s.lutSize - 1);
        check();
        RequireInvalid([&] { Evaluate(evaluator, input, {}); });
        lut[0] = -1;
        RequireInvalid([&] { Evaluate(evaluator, input, lut); });
        lut[0] = s.lutSize;
        RequireInvalid([&] { Evaluate(evaluator, input, lut); });
    }
    Release(s);
    std::cout << " passed\n";
}
} // namespace

int main() {
    try {
        TestSchedules();
        TestCleaning();
        for (uint32_t digits = 1; digits <= 32; ++digits)
            TestTree(digits);
        std::cout << "Trees with 1..32 digits passed\n";
        TestLookup(1, 8, 0, true);
        TestLookup(8, 8, 0, true);
        for (uint32_t bits = 1; bits <= 8; ++bits)
            TestLookup(2 * bits, bits);
        for (uint32_t bits = 3; bits <= 8; ++bits)
            TestLookup(bits, 1);
        TestLookup(5, 2);
        TestLookup(7, 3, 1);
        TestLookup(10, 4, 2);
        TestLookup(13, 4, 3);
        TestLookup(17, 1);
    } catch (const std::exception& error) {
        std::cerr << "\n" << error.what() << '\n';
        return 1;
    }
    std::cout << "Cleartext LUT tests passed\n";
}
