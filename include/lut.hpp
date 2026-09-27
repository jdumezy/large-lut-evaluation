// SPDX-License-Identifier: BSD-2-Clause
#pragma once
#include "digit-decompose.hpp"
#include <span>

namespace large_lut {
// Apply setup.cleaning iterations of 3x^2-2x^3 to the independent selectors,
// then append 1-sum for each digit. Input ciphertexts are not modified.
SplitSelectors PrepareSelectors(const Setup& setup, SplitSelectors selectors);

// Cleartext entries are integers in [0, 2^inputBits), shared by all slots.
// Exactly 2^inputBits entries are required. Input is at inputQ; output is at q.
RLWECiphertext Evaluate(const Setup& setup, const RLWECiphertext& input, std::span<const int64_t> lut);
} // namespace large_lut
