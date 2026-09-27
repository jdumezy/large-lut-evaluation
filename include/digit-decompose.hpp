// SPDX-License-Identifier: BSD-2-Clause
#pragma once
#include "setup.hpp"
#include <span>

namespace large_lut {
using SelectorMatrix = std::vector<std::vector<Ciphertext>>;

struct SplitSelectors {
    SelectorMatrix real;
    SelectorMatrix imag;
};

// Rows run from least to most significant digit. A width-w row contains
// selectors for 0 through 2^w-2; the last selector is omitted. Both matrices
// contain real-valued CKKS slots, covering the two halves of the coefficient
// packing, with bit-reversed order within each half. Recompose followed by
// homomorphic decoding restores coefficient order. Selectors are uncleaned;
// clean before forming the final complement.
// Evaluation uses public/evaluation keys and leaves the input unchanged.
// Optional per-digit level drops are applied to shared powers before polynomial
// evaluation. Empty means no drops. Each drop must fit the computation depth.
SplitSelectors Decompose(const Setup& setup, const RLWECiphertext& input,
                         std::span<const uint32_t> levelDrops = {});

// Join the two real-valued slot vectors as real + i*imag for homomorphic decoding.
Ciphertext Recompose(const Setup& setup, const Ciphertext& real, const Ciphertext& imag);
} // namespace large_lut
