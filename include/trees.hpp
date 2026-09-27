// SPDX-License-Identifier: BSD-2-Clause
#pragma once
#include "digit-decompose.hpp"
#include <span>

namespace large_lut {
// Initial level reductions for the collapsed tree, for 1..32 digit rows.
std::vector<uint32_t> TreeLevelDrops(uint32_t digits);

// Complete selector rows, least significant first, at a common level and scale.
// LUT indexing uses mixed radix: row 0 varies fastest. Its size must equal the
// product of the row widths. Input ciphertexts are not modified.
Ciphertext EvaluateCleartextTree(const Setup& setup, SelectorMatrix selectors, std::span<const int64_t> lut);
} // namespace large_lut
