// SPDX-License-Identifier: BSD-2-Clause
#include "schemelet/rlwe-mp.h"
#include "setup.hpp"
#include <algorithm>
#include <bit>
#include <cmath>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <string>

using namespace large_lut;
using namespace lbcrypto;

namespace {
void Require(bool condition, const std::string& message) {
    if (!condition)
        throw std::runtime_error(message);
}

template <class Function> void RequireInvalid(Function&& function, const std::string& message) {
    try {
        function();
    } catch (const std::invalid_argument&) {
        return;
    }
    throw std::runtime_error("Expected invalid_argument: " + message);
}

void TestLimits() {
    Parameters p;
    for (auto bits : {0u, 33u}) {
        p.inputBits = bits;
        RequireInvalid([&] { ValidateParameters(p); }, "input width");
    }
    p = {};
    for (auto bits : {0u, 9u}) {
        p.digitBits = bits;
        RequireInvalid([&] { ValidateParameters(p); }, "digit width");
    }
    p = {};
    for (auto ring : {128u, 257u}) {
        p.ringDim = ring;
        RequireInvalid([&] { ValidateParameters(p); }, "ring dimension");
    }
    p = {};
    p.toy = false;
    RequireInvalid([&] { ValidateParameters(p); }, "toy ring used in secure mode");
    Require(SecurityBound(65536, 192) == 1481, "Wrong security bound for h=192");
    RequireInvalid([] { SecurityBound(131072, 192); }, "extrapolated security bound");
    p = {};
    p.levelBudget = {0, 2};
    RequireInvalid([&] { ValidateParameters(p); }, "empty encoding budget");
    p.levelBudget = {2};
    RequireInvalid([&] { ValidateParameters(p); }, "incomplete budget");
    p = {};
    p.inputBits = 32;
    p.digitBits = 8;
    ValidateParameters(p);
    Require(DigitWidths(32, 8) == std::vector<uint32_t>({8, 8, 8, 8}), "32-bit domain");
    Require(DigitWidths(10, 4) == std::vector<uint32_t>({4, 4, 2}), "partial top digit");
    Require(DigitWidths(32, 7) == std::vector<uint32_t>({7, 7, 7, 7, 4}), "32-bit partial digit");
    Require(DigitWidths(1, 8) == std::vector<uint32_t>({1}), "digit wider than domain");
    Require(SecurityBound(4096, 32) == 19, "Security table lower corner");
    Require(SecurityBound(65536, 1024) == 1714, "Security table upper corner");
    Require(SecurityBound(65536, 256) == 1553, "Security bound for h=256");
    RequireInvalid([] { SecurityBound(65536, 200); }, "unlisted weight");
    p = {};
    p.h = 0;
    RequireInvalid([&] { ValidateParameters(p); }, "zero weight");
    p.h = p.ringDim + 1;
    RequireInvalid([&] { ValidateParameters(p); }, "weight larger than ring");

    p = {};
    Require(SelectLevelBudget(p) == std::vector<uint32_t>({4, 3}), "8-bit LUT budget");
    p.inputBits = 16;
    Require(SelectLevelBudget(p) == std::vector<uint32_t>({3, 2}), "16-bit LUT budget");
    p.inputBits = 32;
    Require(SelectLevelBudget(p) == std::vector<uint32_t>({4, 2}), "32-bit LUT budget");
    p.levelBudget = {2, 1};
    Require(SelectLevelBudget(p) == p.levelBudget, "Explicit level budget override");
}

void TestSecurityRejection() {
    Parameters p;
    p.inputBits = 1;
    p.digitBits = 1;
    p.ringDim = 32768;
    p.levelBudget = {2, 2};
    p.toy = false;
    try {
        CreateSetup(p);
    } catch (const std::invalid_argument& error) {
        Require(std::string(error.what()).starts_with("log2(QP)="),
                "Expected rejection by the QP security check");
        CryptoContextFactory<DCRTPoly>::ReleaseAllContexts();
        return;
    }
    throw std::runtime_error("Insecure QP configuration accepted");
}

double HermiteValue(const Coefficients& coefficients, uint32_t x, uint32_t base) {
    const auto z = std::polar(1.0, 2.0 * std::numbers::pi * x / base);
    std::complex<double> value{};
    for (auto it = coefficients.rbegin(); it != coefficients.rend(); ++it)
        value = value * z + *it;
    // Coefficients already include the 1/2 needed by f(z) + conjugate(f(z)).
    return 2.0 * value.real();
}

void TestCoefficients(const Setup& s) {
    const uint32_t base = 1u << s.digitBits;
    Require(s.selectorCoefficients.size() == base - 1, "Omitted selector count");
    for (uint32_t x = 0; x < base; ++x) {
        const double digit = s.digitBits == 1
                                 ? (s.modCoefficients[0] + s.modCoefficients[1] * double(1 - x)).real()
                                 : HermiteValue(s.modCoefficients, x, base) * s.modScale;
        Require(std::abs(digit - x) < 1e-8, "Digit coefficient scaling");
        double sum = 0;
        for (uint32_t j = 0; j < base - 1; ++j) {
            const double selector =
                s.digitBits == 1
                    ? (s.selectorCoefficients[j][0] + s.selectorCoefficients[j][1] * double(1 - x)).real()
                    : HermiteValue(s.selectorCoefficients[j], x, base);
            Require(std::abs(selector - double(x == j)) < 1e-9, "Selector coefficient scaling");
            sum += selector;
        }
        Require(std::abs(1.0 - sum - double(x == base - 1)) < 1e-8, "Complement selector");
    }
    // A shorter final digit is embedded at multiples of base/topBase.
    const uint32_t topBase = 1u << s.digitWidths.back();
    if (topBase != base) {
        for (uint32_t x = 0; x < topBase; ++x) {
            double sum = 0;
            for (uint32_t j = 0; j < topBase - 1; ++j) {
                const double value =
                    HermiteValue(s.selectorCoefficients[j * (base / topBase)], x * (base / topBase), base);
                Require(std::abs(value - double(x == j)) < 1e-9, "Partial digit coefficient mapping");
                sum += value;
            }
            Require(std::abs(1.0 - sum - double(x == topBase - 1)) < 1e-8, "Partial complement");
        }
    }
}

void TestToy(Parameters p) {
    auto s = CreateSetup(p);
    Require(s.parameters.toy, "Tests must use toy parameters");
    Require(s.lutSize == (uint64_t{1} << p.inputBits), "LUT count overflow");
    Require(s.inputQ * (BigInteger(1) << s.digitBits) == s.q * s.inputModulus,
            "Initial extraction scale invariant");
    TestCoefficients(s);

    const auto crypto = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(s.cc->GetCryptoParameters());
    double logQ = 0, logP = 0;
    for (const auto& prime : crypto->GetElementParams()->GetParams())
        logQ += std::log2(prime->GetModulus().ConvertToDouble());
    for (const auto& prime : crypto->GetParamsP()->GetParams())
        logP += std::log2(prime->GetModulus().ConvertToDouble());
    Require(logP > 0 && std::abs(s.logQP - logQ - logP) < 1e-8, "QP must include auxiliary P");
    auto secret = s.secretKey->GetPrivateElement().GetElementAtIndex(0);
    secret.SetFormat(Format::COEFFICIENT);
    size_t weight = 0;
    for (size_t i = 0; i < secret.GetLength(); ++i)
        weight += secret[i] != NativeInteger(0);
    Require(weight == p.h, "Incorrect secret Hamming weight");

    const std::vector<double> slots{0.25, -0.5, 1.0};
    auto plaintext = s.cc->MakeCKKSPackedPlaintext(slots);
    auto encrypted = s.cc->Encrypt(s.publicKey, plaintext);
    Plaintext decrypted;
    s.cc->Decrypt(s.secretKey, encrypted, &decrypted);
    const auto decoded = decrypted->GetRealPackedValue();
    for (size_t i = 0; i < slots.size(); ++i)
        Require(std::abs(decoded[i] - slots[i]) < 1e-8, "Public/secret key mismatch");

    std::vector<int64_t> values(p.ringDim);
    for (size_t i = 0; i < values.size(); ++i)
        values[i] = (uint64_t(i) * 2654435761ULL) % s.lutSize;
    values[0] = 0;
    values[1] = s.lutSize - 1;
    values[values.size() / 2] = s.lutSize / 2;
    values.back() = s.lutSize - 1;
    auto ct = EncryptInput(s, values);
    Require(DecryptInput(s, ct) == values, "Full packed RLWE input round trip");
    // The built-in FBT exercises context precomputation, keys, coefficient
    // normalization, and decoding without depending on iterative extraction.
    auto digitInput = ct;
    for (auto& poly : digitInput)
        poly.SwitchModulus(s.q, 1, 0, 0);
    auto ckks = SchemeletRLWEMP::ConvertRLWEToCKKS(*s.cc, digitInput, s.publicKey, s.q, s.slots, s.depth);
    auto digit =
        s.cc->EvalFBT(ckks, s.modCoefficients, s.digitBits, s.qPrime, s.modScale, s.computationDepth);
    auto digitRLWE = SchemeletRLWEMP::ConvertCKKSToRLWE(digit, s.q);
    auto expectedDigits = values;
    for (auto& value : expectedDigits)
        value &= (uint64_t{1} << s.digitBits) - 1;
    Require(DecryptOutput(s, digitRLWE) == expectedDigits, "FBT digit scaling round trip");

    // Exercise output decoding independently, at the output modulus q.
    auto output = SchemeletRLWEMP::EncryptCoeff(values, s.q, s.outputModulus, s.secretKey, s.elementParams);
    Require(DecryptOutput(s, output) == values, "Full packed RLWE output round trip");
    RequireInvalid([&] { EncryptInput(s, {0, 1}); }, "short input vector");
    values[0] = -1;
    RequireInvalid([&] { EncryptInput(s, values); }, "negative input");
    values[0] = s.lutSize;
    RequireInvalid([&] { EncryptInput(s, values); }, "out-of-domain input");
    RequireInvalid([&] { DecryptInput(s, {}); }, "malformed ciphertext");
    if (s.inputQ != s.q)
        RequireInvalid([&] { DecryptOutput(s, ct); }, "input ciphertext at output modulus");
    std::cout << "toy: inputBits=" << p.inputBits << " digitBits=" << p.digitBits << " h=" << p.h
              << " cleaning=" << s.cleaning << " depth=" << s.depth << " logQP=" << s.logQP << '\n';
    s.cc->ClearEvalMultKeys();
    s.cc->ClearEvalAutomorphismKeys();
    CryptoContextFactory<DCRTPoly>::ReleaseAllContexts();
}
} // namespace

int main() {
    try {
        TestLimits();
        TestSecurityRejection();
        // Each supported digit size, with a partial final digit where possible.
        for (uint32_t bits = 1; bits <= 8; ++bits) {
            Parameters p;
            p.inputBits = 2 * bits + 1;
            p.digitBits = bits;
            TestToy(p);
        }
        Parameters p;
        p.inputBits = 32;
        p.digitBits = 8;
        TestToy(p); // No LUT allocation or large-domain evaluation.
        p = {};
        p.inputBits = 1;
        p.digitBits = 8;
        p.lutKind = LutKind::Ciphertext;
        TestToy(p);
        for (auto h : {32u, 64u, 128u, 256u, 512u, 1024u}) {
            p = {};
            p.h = h;
            p.ringDim = std::max(256u, std::bit_ceil(h));
            TestToy(p);
        }
        p = {};
        p.lutKind = LutKind::Plaintext;
        TestToy(p);
        std::cout << "setup tests passed\n";
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
