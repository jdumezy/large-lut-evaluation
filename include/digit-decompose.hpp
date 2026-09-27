// SPDX-License-Identifier: BSD-2-Clause
#pragma once
#include "setup.hpp"

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
SplitSelectors Decompose(const Setup& setup, const RLWECiphertext& input);

// Join the two real-valued slot vectors as real + i*imag for homomorphic decoding.
Ciphertext Recompose(const Setup& setup, const Ciphertext& real, const Ciphertext& imag);
} // namespace large_lut
