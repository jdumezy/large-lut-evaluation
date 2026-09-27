// SPDX-License-Identifier: BSD-2-Clause
#pragma once
#include "openfhe.h"
#include <complex>
#include <cstdint>
#include <vector>

namespace large_lut {
using Ciphertext = lbcrypto::Ciphertext<lbcrypto::DCRTPoly>;
using RLWECiphertext = std::vector<lbcrypto::Poly>;
using Coefficients = std::vector<std::complex<double>>;

enum class LutKind {
    Cleartext, // Integer constants shared by all slots.
    Plaintext, // Encoded CKKS plaintexts, potentially different in each slot.
    Ciphertext // Encrypted CKKS entries.
};
struct Parameters {
    uint32_t inputBits = 8;
    uint32_t digitBits = 4;
    uint32_t ringDim = 256; // Toy default; secure mode requires a table-covered N.
    bool toy = true;
    LutKind lutKind = LutKind::Cleartext;
    uint32_t h = 192;
    int cleaning = -1;      // -1: 0 up to 12 bits, 1 up to 20, 2 at 21, 3 above.
    uint32_t scaleBits = 0; // 0: 48 up to 18 input bits, 59 above.
    uint32_t largeDigits = 3;
    std::vector<uint32_t> levelBudget; // Empty: select from inputBits and digitBits.
};

struct Setup {
    Parameters parameters;
    uint64_t lutSize = 0;
    std::vector<uint32_t> digitWidths;
    uint32_t digitBits = 0, slots = 0, depth = 0, computationDepth = 0;
    uint32_t cleaning = 0, modScale = 1;
    double logQP = 0;
    lbcrypto::BigInteger inputModulus, outputModulus, q, inputQ, encryptionQ, qPrime;
    Coefficients modCoefficients;
    std::vector<Coefficients> selectorCoefficients; // Last full-base selector omitted.
    lbcrypto::CryptoContext<lbcrypto::DCRTPoly> cc;
    lbcrypto::PublicKey<lbcrypto::DCRTPoly> publicKey;
    lbcrypto::PrivateKey<lbcrypto::DCRTPoly> secretKey;
    std::shared_ptr<lbcrypto::ILDCRTParams<lbcrypto::DCRTPoly::Integer>> elementParams;
};

// Preflight validation does not create a crypto context or allocate a LUT.
void ValidateParameters(const Parameters& parameters);
std::vector<uint32_t> DigitWidths(uint32_t inputBits, uint32_t digitBits);
std::vector<uint32_t> SelectLevelBudget(const Parameters& parameters);
double SecurityBound(uint32_t ringDim, uint32_t h);
Setup CreateSetup(Parameters parameters = {});
// Inputs are coefficient-packed integers in [0, 2^inputBits); exactly ringDim
// values.
RLWECiphertext EncryptInput(const Setup& setup, const std::vector<int64_t>& values);
std::vector<int64_t> DecryptInput(const Setup& setup, const RLWECiphertext& ciphertext);
std::vector<int64_t> DecryptOutput(const Setup& setup, const RLWECiphertext& ciphertext);
} // namespace large_lut
