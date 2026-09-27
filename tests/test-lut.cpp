// SPDX-License-Identifier: BSD-2-Clause
#include "lut.hpp"
#include "trees.hpp"
#include <algorithm>
#include <bit>
#include <cmath>
#include <iostream>
#include <limits>
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
    const auto drops = TreeLevelDrops(digits);
    for (size_t digit = 0; digit < selectors.size(); ++digit)
        for (auto& ct : selectors[digit])
            if (drops[digit])
                ct = s.cc->LevelReduce(ct, nullptr, drops[digit]);
    auto earlyResult = EvaluateCleartextTree(s, selectors, lut);
    auto earlyValues = Slots(s, earlyResult);
    Require(earlyResult->GetLevel() == result->GetLevel(), "Early drops changed the tree depth");
    for (size_t slot = 0; slot < s.slots; ++slot)
        Require(std::abs(earlyValues[slot] - double(lut[indices[slot]])) < 1e-6,
                "Tree with early drops lookup mismatch");
    selectors.back()[0] = s.cc->LevelReduce(selectors.back()[0], nullptr, 1);
    RequireInvalid([&] { EvaluateCleartextTree(s, selectors, lut); });
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
        RequireInvalid([&] { Evaluate(evaluator, input, std::span<const int64_t>{}); });
        lut[0] = -1;
        RequireInvalid([&] { Evaluate(evaluator, input, lut); });
        lut[0] = s.lutSize;
        RequireInvalid([&] { Evaluate(evaluator, input, lut); });
    }
    Release(s);
    std::cout << " passed\n";
}
void TestIntegerContraction() {
    auto p = ToyParameters(8, 4, 0);
    p.scaleBits = 59;
    auto s = CreateSetup(p);
    std::mt19937_64 random(4061);
    std::vector<Ciphertext> row;
    for (size_t i = 0; i < 256; ++i) {
        std::vector<double> values(s.slots);
        for (auto& value : values)
            value = double(random() % 1000) / 1000;
        row.push_back(s.cc->Encrypt(s.publicKey, s.cc->MakeCKKSPackedPlaintext(values)));
    }
    const auto unchanged = row[0]->GetElements();
    for (uint32_t size : {1u, 2u, 16u, 256u}) {
        std::vector<Ciphertext> selected(row.begin(), row.begin() + size);
        for (uint64_t maximum : {uint64_t{0}, uint64_t{65535}, uint64_t{0xffffffff},
                                 uint64_t(std::numeric_limits<int64_t>::max())}) {
            std::vector<int64_t> lut(size);
            for (auto& value : lut)
                value = maximum ? random() % maximum : 0;
            lut[0] = maximum;
            auto expected = s.cc->GetScheme()->MultByInteger(selected[0], lut[0]);
            for (size_t i = 1; i < lut.size(); ++i)
                s.cc->EvalAddInPlaceNoCheck(expected, s.cc->GetScheme()->MultByInteger(selected[i], lut[i]));
            auto actual = EvaluateCleartextTree(s, {selected}, lut);
            Require(actual->GetElements() == expected->GetElements(),
                    "Integer contraction changed ciphertext coefficients");
            Require(actual->GetLevel() == expected->GetLevel() &&
                        actual->GetNoiseScaleDeg() == expected->GetNoiseScaleDeg() &&
                        actual->GetScalingFactor() == expected->GetScalingFactor(),
                    "Integer contraction changed scale metadata");
        }
    }
    Require(row[0]->GetElements() == unchanged, "Integer contraction mutated a selector");
    Release(s);
    std::cout << "Exact integer contraction checks passed\n";
}

void TestEncodedLookup(uint32_t bits, uint32_t digitBits, int cleaning, uint32_t scaleBits = 48) {
    std::cout << "Encoded LUT " << bits << "/" << digitBits << ", cleaning " << cleaning << std::flush;
    std::mt19937_64 random(903 + bits);
    std::vector<int64_t> table(uint64_t{1} << bits);
    for (auto& value : table)
        value = random() % table.size();
    for (auto kind : {LutKind::Cleartext, LutKind::Plaintext, LutKind::Ciphertext}) {
        auto p = ToyParameters(bits, digitBits, cleaning);
        p.lutKind = kind;
        p.scaleBits = scaleBits;
        auto s = CreateSetup(p);
        auto evaluator = s;
        evaluator.secretKey.reset();
        std::vector<int64_t> values(p.ringDim);
        for (size_t i = 0; i < values.size(); ++i)
            values[i] = i % s.lutSize;
        auto verify = [&](const auto& lut, bool varying) {
            for (size_t repetition = 0; repetition < 2; ++repetition) {
                std::shuffle(values.begin(), values.end(), random);
                auto input = EncryptInput(s, values);
                const auto saved = input;
                auto actual = DecryptOutput(s, Evaluate(evaluator, input, lut));
                Require(input == saved, "Encoded LUT evaluation changed the RLWE input");
                for (size_t i = 0; i < values.size(); ++i) {
                    auto expected = table[values[i]];
                    if (varying)
                        expected = (expected + 7 * i + 31 * (i / s.slots)) % s.lutSize;
                    Require(actual[i] == expected,
                            "Encoded LUT mismatch at coefficient " + std::to_string(i));
                }
            }
        };
        if (kind == LutKind::Cleartext) {
            verify(table, false);
            RequireInvalid([&] { EncodeLut(s, table); });
        } else {
            for (bool varying : {false, true}) {
                auto plain = varying ? EncodeLut(s,
                                                 [&](uint64_t index, uint32_t coefficient) -> int64_t {
                                                     return (table[index] + 7 * coefficient +
                                                             31 * (coefficient / s.slots)) %
                                                            s.lutSize;
                                                 })
                                     : EncodeLut(s, table);
                auto savedReal = plain.real[0]->GetElement<DCRTPoly>();
                auto savedImag = plain.imag.back()->GetElement<DCRTPoly>();
                if (kind == LutKind::Plaintext) {
                    verify(plain, varying);
                    auto malformed = plain;
                    malformed.imag.pop_back();
                    RequireInvalid([&] { Evaluate(evaluator, {}, malformed); });
                    malformed = plain;
                    malformed.real[0] = nullptr;
                    RequireInvalid([&] { Evaluate(evaluator, {}, malformed); });
                    auto level = plain.real[0]->GetLevel();
                    plain.real[0]->SetLevel(level + 1);
                    RequireInvalid([&] { Evaluate(evaluator, {}, plain); });
                    plain.real[0]->SetLevel(level);
                    RequireInvalid([&] { EncryptLut(evaluator, plain); });
                } else {
                    auto encrypted = EncryptLut(evaluator, plain);
                    auto savedCipherReal = encrypted.real[0]->GetElements();
                    auto savedCipherImag = encrypted.imag.back()->GetElements();
                    verify(encrypted, varying);
                    Require(encrypted.real[0]->GetElements() == savedCipherReal &&
                                encrypted.imag.back()->GetElements() == savedCipherImag,
                            "Evaluation mutated encrypted LUT entries");
                    auto malformed = encrypted;
                    malformed.real[0] = malformed.real[0]->Clone();
                    malformed.real[0]->SetKeyTag("different-key");
                    RequireInvalid([&] { Evaluate(evaluator, {}, malformed); });
                    malformed = encrypted;
                    malformed.imag[0] = s.cc->LevelReduce(malformed.imag[0], nullptr, 1);
                    RequireInvalid([&] { Evaluate(evaluator, {}, malformed); });
                    RequireInvalid([&] { Evaluate(evaluator, {}, plain); });
                }
                Require(plain.real[0]->GetElement<DCRTPoly>() == savedReal &&
                            plain.imag.back()->GetElement<DCRTPoly>() == savedImag,
                        "Evaluation or encryption mutated plaintext LUT entries");
            }
            RequireInvalid([&] { EncodeLut(s, [](uint64_t, uint32_t) { return -1; }); });
            RequireInvalid([&] { EncodeLut(s, [&](uint64_t, uint32_t) { return int64_t(s.lutSize); }); });
            RequireInvalid([&] { Evaluate(evaluator, {}, std::span<const int64_t>(table)); });
        }
        Release(s);
    }
    std::cout << " passed\n";
}
} // namespace

int main() {
    try {
        TestIntegerContraction();
        TestEncodedLookup(1, 8, 0);
        TestEncodedLookup(4, 4, 0);
        TestEncodedLookup(4, 2, 1);
        TestEncodedLookup(5, 2, 0);
        TestEncodedLookup(7, 2, 2);
        TestEncodedLookup(7, 3, 1);
        TestEncodedLookup(8, 2, 0);
        TestEncodedLookup(5, 1, 3, 59);
        TestEncodedLookup(9, 2, 1, 59);
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
        TestLookup(15, 4);
        TestLookup(16, 4);
        TestLookup(17, 4);
        TestLookup(17, 1);
    } catch (const std::exception& error) {
        std::cerr << "\n" << error.what() << '\n';
        return 1;
    }
    std::cout << "Cleartext LUT tests passed\n";
}
