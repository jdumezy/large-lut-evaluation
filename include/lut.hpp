// SPDX-License-Identifier: BSD-2-Clause
#pragma once
#include "digit-decompose.hpp"
#include <functional>
#include <span>

namespace large_lut {
// Apply setup.cleaning iterations of 3x^2-2x^3 to the independent selectors,
// then append 1-sum for each digit. Input ciphertexts are not modified.
SplitSelectors PrepareSelectors(const Setup& setup, SplitSelectors selectors);

// Cleartext entries are integers in [0, 2^inputBits), shared by all slots.
// Exactly 2^inputBits entries are required. Input is at inputQ; output is at q.
RLWECiphertext Evaluate(const Setup& setup, const RLWECiphertext& input, std::span<const int64_t> lut);
struct PlaintextLut {
    std::vector<lbcrypto::Plaintext> real, imag;
};
struct CiphertextLut {
    std::vector<Ciphertext> real, imag;
};

// The callback returns the entry for (LUT index, input coefficient index).
// Values are in [0, 2^inputBits). Packing and entry levels are selected internally.
using LutValues = std::function<int64_t(uint64_t, uint32_t)>;
PlaintextLut EncodeLut(const Setup& setup, const LutValues& values);
PlaintextLut EncodeLut(const Setup& setup, std::span<const int64_t> values);
CiphertextLut EncryptLut(const Setup& setup, const PlaintextLut& lut);

RLWECiphertext Evaluate(const Setup& setup, const RLWECiphertext& input, const PlaintextLut& lut);
RLWECiphertext Evaluate(const Setup& setup, const RLWECiphertext& input, const CiphertextLut& lut);
} // namespace large_lut
