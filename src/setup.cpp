// SPDX-License-Identifier: BSD-2-Clause
#include "setup.hpp"
#include "math/hermite.h"
#include "schemelet/rlwe-mp.h"
#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <stdexcept>

namespace large_lut {
using namespace lbcrypto;

void ValidateParameters(const Parameters& p) {
    if (p.inputBits < 1 || p.inputBits > 32 || p.digitBits < 1 || p.digitBits > 8)
        throw std::invalid_argument("inputBits must be 1..32 and digitBits 1..8");
    if (!std::has_single_bit(p.ringDim) || p.ringDim < 256)
        throw std::invalid_argument("ringDim must be a power of two >= 256");
    if (p.h < 1 || p.h > p.ringDim)
        throw std::invalid_argument("h must be in 1..ringDim");
    if (p.cleaning < -1 || p.cleaning > 8)
        throw std::invalid_argument("cleaning must be -1 (automatic) or 0..8");
    if (p.scaleBits && (p.scaleBits < 40 || p.scaleBits > 59))
        throw std::invalid_argument("scaleBits must be 0 (automatic) or 40..59");
    if (p.lutKind != LutKind::Cleartext && p.lutKind != LutKind::Plaintext &&
        p.lutKind != LutKind::Ciphertext)
        throw std::invalid_argument("Invalid LUT kind");
    if (p.largeDigits < 1)
        throw std::invalid_argument("largeDigits must be positive");
    const uint32_t logSlots = std::bit_width(p.ringDim / 2) - 1;
    if (!p.levelBudget.empty() &&
        (p.levelBudget.size() != 2 || p.levelBudget[0] < 1 || p.levelBudget[1] < 1 ||
         p.levelBudget[0] > logSlots || p.levelBudget[1] > logSlots))
        throw std::invalid_argument("levelBudget needs two values in 1..log2(slots)");
    if (!p.toy) {
        SecurityBound(p.ringDim, p.h);
        SecurityBound(p.ringDim, 32);
    }
}

std::vector<uint32_t> DigitWidths(uint32_t inputBits, uint32_t digitBits) {
    if (inputBits < 1 || inputBits > 32 || digitBits < 1 || digitBits > 8)
        throw std::invalid_argument("inputBits must be 1..32 and digitBits 1..8");
    std::vector<uint32_t> widths;
    for (uint32_t offset = 0; offset < inputBits; offset += digitBits)
        widths.push_back(std::min(digitBits, inputBits - offset));
    return widths;
}

std::vector<uint32_t> SelectLevelBudget(const Parameters& p) {
    ValidateParameters(p);
    if (!p.levelBudget.empty())
        return p.levelBudget;

    const auto widths = DigitWidths(p.inputBits, p.digitBits);
    const auto digits = widths.size();
    const auto bits = widths.front();
    std::vector<uint32_t> budget{4, 2};
    // Defaults balance transform cost against the levels consumed by
    // interpolation, cleaning, and the multiplexer tree.
    switch (bits) {
    case 1:
        budget = digits == 2 ? std::vector<uint32_t>{4, 4} : std::vector<uint32_t>{3, 2};
        break;
    case 2:
        if (digits <= 4)
            budget = {4, 4};
        else if (digits == 5)
            budget = {4, 3};
        else if (digits == 6)
            budget = {3, 3};
        else if (digits == 7)
            budget = {3, 2};
        else if (digits <= 9)
            budget = {2, 2};
        break;
    case 3:
        if (digits <= 3)
            budget = {4, 4};
        else if (digits == 4)
            budget = {4, 3};
        else if (digits == 5)
            budget = {3, 2};
        else if (digits == 6)
            budget = {2, 2};
        break;
    case 4:
        if (digits <= 3)
            budget = {4, 3};
        else if (p.inputBits <= 18)
            budget = {3, 2};
        break;
    case 5:
        if (digits <= 2)
            budget = {4, 3};
        else if (digits == 3)
            budget = p.inputBits <= 12 ? std::vector<uint32_t>{3, 3} : std::vector<uint32_t>{2, 2};
        else if (p.inputBits <= 18)
            budget = {3, 2};
        break;
    case 6:
        if (digits <= 2)
            budget = {3, 3};
        else if (digits == 3)
            budget = {2, 2};
        break;
    case 7:
        if (digits <= 2)
            budget = {3, 2};
        else if (digits == 3)
            budget = {2, 2};
        break;
    case 8:
        if (digits <= 2)
            budget = p.inputBits < 12 ? std::vector<uint32_t>{3, 3} : std::vector<uint32_t>{2, 2};
        break;
    }
    const uint32_t logSlots = std::bit_width(p.ringDim / 2) - 1;
    for (auto& levels : budget)
        levels = std::min(levels, logSlots);
    return budget;
}

double SecurityBound(uint32_t ringDim, uint32_t h) {
    // https://github.com/jdumezy/sparse-key-estimate/blob/master/Precomputed-Tables/128bits_security.md
    constexpr std::array<uint32_t, 7> weights{32, 64, 128, 192, 256, 512, 1024};
    constexpr uint32_t bounds[5][7] = {
        {19, 58, 86, 96, 101, 106, 111},          // log2(N) = 12
        {37, 113, 167, 187, 197, 209, 214},       // log2(N) = 13
        {77, 229, 333, 370, 388, 416, 426},       // log2(N) = 14
        {164, 468, 669, 739, 774, 829, 854},      // log2(N) = 15
        {349, 961, 1351, 1481, 1553, 1659, 1714}, // log2(N) = 16
    };
    const auto column = std::find(weights.begin(), weights.end(), h);
    if (!std::has_single_bit(ringDim) || ringDim < 4096 || ringDim > 65536 || column == weights.end())
        throw std::invalid_argument("No tabulated 128-bit security bound for this (ringDim, h)");
    return bounds[std::bit_width(ringDim) - 13][column - weights.begin()];
}

namespace {
KeyPair<DCRTPoly> GenerateKeys(const CryptoContext<DCRTPoly>& cc, uint32_t h) {
    if (h == 192)
        return cc->KeyGen();

    const auto crypto = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(cc->GetCryptoParameters());
    const auto params = crypto->GetParamsPK();
    DCRTPoly::TugType ternary;
    DCRTPoly::DugType uniform;
    DCRTPoly secret(ternary, params, Format::EVALUATION, h);
    DCRTPoly a(uniform, params, Format::EVALUATION);
    DCRTPoly error(crypto->GetDiscreteGaussianGenerator(), params, Format::EVALUATION);
    auto b = error * crypto->GetNoiseScale() - a * secret;

    const auto extraTowers = params->GetParams().size() - crypto->GetElementParams()->GetParams().size();
    if (extraTowers)
        secret.DropLastElements(extraTowers);
    KeyPair<DCRTPoly> keys(std::make_shared<PublicKeyImpl<DCRTPoly>>(cc),
                           std::make_shared<PrivateKeyImpl<DCRTPoly>>(cc));
    keys.secretKey->SetPrivateElement(std::move(secret));
    keys.publicKey->SetPublicElements({std::move(b), std::move(a)});
    keys.publicKey->SetKeyTag(keys.secretKey->GetKeyTag());
    return keys;
}
} // namespace

Setup CreateSetup(Parameters p) {
    ValidateParameters(p);
    p.levelBudget = SelectLevelBudget(p);
    Setup s;
    s.parameters = p;
    s.lutSize = uint64_t{1} << p.inputBits;
    s.digitWidths = DigitWidths(p.inputBits, p.digitBits);
    s.digitBits = s.digitWidths.front();
    const uint32_t base = uint32_t{1} << s.digitBits;
    s.slots = p.ringDim / 2;
    s.cleaning = p.cleaning >= 0 ? p.cleaning
                                 : (p.inputBits <= 12   ? 0
                                    : p.inputBits <= 20 ? 1
                                    : p.inputBits == 21 ? 2
                                                        : 3);
    s.parameters.cleaning = s.cleaning;
    s.parameters.scaleBits = p.scaleBits ? p.scaleBits : (p.inputBits <= 18 ? 48 : 59);
    s.inputModulus = BigInteger(1) << p.inputBits;
    s.outputModulus = s.inputModulus;
    s.q = BigInteger(1) << s.parameters.scaleBits;
    // Invariant at the first extraction: inputQ / P = q / base.
    s.inputQ = s.q << (p.inputBits - s.digitBits);
    s.encryptionQ = BigInteger(1) << (s.parameters.scaleBits + p.inputBits + 16);
    s.modScale = s.digitBits > 1 ? uint32_t{1} << (s.digitBits - 2) : 1;
    if (s.digitBits == 1) {
        // OpenFHE's binary FBT shortcut uses an affine polynomial in its
        // precomputed cosine square, rather than the general Hermite coefficients.
        s.modCoefficients = {1., -1.};
        s.selectorCoefficients = {{0., 1.}};
    } else {
        s.modCoefficients =
            GetHermiteTrigCoefficients([base](int64_t x) { return x % base; }, base, 1, s.modScale);
        for (uint32_t j = 0; j < base - 1; ++j)
            s.selectorCoefficients.push_back(
                GetHermiteTrigCoefficients([j](int64_t x) -> int64_t { return x == j; }, base, 1, 1.));
    }

    // A balanced cleartext tree needs ceil(log2(number of digits)) levels.
    // Encoded plaintext and ciphertext entries reserve one extra multiplication level. Each cleaning
    // pass (3x^2 - 2x^3) consumes two. Decoding is included in GetFBTDepth.
    const auto treeDepth = std::bit_width(s.digitWidths.size() - 1);
    s.computationDepth = treeDepth + (p.lutKind != LutKind::Cleartext) + 2 * s.cleaning;
    // Modulus raising uses a separate weight-32 key, keeping the approximation
    // interval independent of the main secret's Hamming weight.
    constexpr auto secretDistribution = SPARSE_ENCAPSULATED;
    s.depth = s.computationDepth + FHECKKSRNS::GetFBTDepth(p.levelBudget, s.modCoefficients, BigInteger(base),
                                                           1, secretDistribution);
    if (p.largeDigits > s.depth + 1)
        throw std::invalid_argument("largeDigits exceeds the number of Q towers");
    CCParams<CryptoContextCKKSRNS> cp;
    cp.SetSecretKeyDist(secretDistribution);
    cp.SetSecurityLevel(HEStd_NotSet); // Sparse-key bound checked explicitly below.
    cp.SetScalingTechnique(FIXEDMANUAL);
    cp.SetKeySwitchTechnique(HYBRID);
    cp.SetScalingModSize(s.parameters.scaleBits);
    cp.SetFirstModSize(s.parameters.scaleBits);
    cp.SetNumLargeDigits(p.largeDigits);
    cp.SetRingDim(p.ringDim);
    cp.SetBatchSize(s.slots);
    cp.SetMultiplicativeDepth(s.depth);
    cp.SetMaxRelinSkDeg(3); // Lazy tree products can have degree three.
    s.cc = GenCryptoContext(cp);
    for (auto feature : {PKE, KEYSWITCH, LEVELEDSHE, ADVANCEDSHE, FHE})
        s.cc->Enable(feature);
    const auto crypto = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(s.cc->GetCryptoParameters());
    for (const auto& prime : crypto->GetParamsQP()->GetParams())
        s.logQP += std::log2(prime->GetModulus().ConvertToDouble());
    if (!p.toy) {
        if (s.logQP > SecurityBound(p.ringDim, p.h))
            throw std::invalid_argument("log2(QP)=" + std::to_string(s.logQP) +
                                        " exceeds the 128-bit security bound for h=" + std::to_string(p.h));
        const auto& q0 = crypto->GetElementParams()->GetParams().front()->GetModulus();
        const auto& p0 = crypto->GetParamsP()->GetParams().front()->GetModulus();
        const double logSparseQP = std::log2(q0.ConvertToDouble()) + std::log2(p0.ConvertToDouble());
        if (logSparseQP > SecurityBound(p.ringDim, 32))
            throw std::invalid_argument("Sparse encapsulation modulus exceeds the 128-bit security bound");
    }
    auto keys = GenerateKeys(s.cc, p.h);
    s.publicKey = keys.publicKey;
    s.secretKey = keys.secretKey;
    s.cc->EvalMultKeysGen(s.secretKey);
    s.cc->EvalFBTSetup(s.modCoefficients, s.slots, BigInteger(base), s.outputModulus, s.q, s.publicKey,
                       {0, 0}, p.levelBudget, 0, s.computationDepth, 1);
    s.cc->EvalBootstrapKeyGen(s.secretKey, s.slots);
    s.elementParams = SchemeletRLWEMP::GetElementParams(s.secretKey, s.depth);
    s.qPrime = s.elementParams->GetModulus();
    return s;
}

RLWECiphertext EncryptInput(const Setup& s, const std::vector<int64_t>& values) {
    if (values.size() != s.parameters.ringDim)
        throw std::invalid_argument("EncryptInput requires exactly ringDim integer values");
    for (auto x : values)
        if (x < 0 || uint64_t(x) >= s.lutSize)
            throw std::invalid_argument("Input integer outside the LUT domain");
    auto ct =
        SchemeletRLWEMP::EncryptCoeff(values, s.encryptionQ, s.inputModulus, s.secretKey, s.elementParams);
    // ModSwitch takes the destination modulus before the source modulus.
    SchemeletRLWEMP::ModSwitch(ct, s.inputQ, s.encryptionQ);
    return ct;
}

namespace {
std::vector<int64_t> DecryptIntegers(const Setup& s, const RLWECiphertext& ct, const BigInteger& modulus) {
    if (ct.size() != 2)
        throw std::invalid_argument("An RLWE ciphertext must have two polynomials");
    for (const auto& poly : ct)
        if (poly.GetRingDimension() != s.parameters.ringDim || poly.GetModulus() != modulus)
            throw std::invalid_argument("Unexpected RLWE ring dimension or modulus");
    auto values = SchemeletRLWEMP::DecryptCoeff(ct, modulus, s.inputModulus, s.secretKey, s.elementParams,
                                                s.slots, s.parameters.ringDim);
    for (auto& x : values) {
        x %= static_cast<int64_t>(s.lutSize);
        if (x < 0)
            x += s.lutSize;
    }
    return values;
}
} // namespace

std::vector<int64_t> DecryptInput(const Setup& s, const RLWECiphertext& ct) {
    return DecryptIntegers(s, ct, s.inputQ);
}

std::vector<int64_t> DecryptOutput(const Setup& s, const RLWECiphertext& ct) {
    return DecryptIntegers(s, ct, s.q);
}
} // namespace large_lut
