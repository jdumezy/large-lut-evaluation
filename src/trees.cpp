// SPDX-License-Identifier: BSD-2-Clause
#include "trees.hpp"
#include <bit>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace large_lut {
namespace {
using Row = std::vector<Ciphertext>;

// Each root side is assembled from power-of-two blocks. The last block may
// be shorter; its leaves start at lower levels to meet the longer branch.
void FillDrops(std::span<uint32_t> drops, uint32_t depth) {
    if (drops.size() == 1) {
        drops[0] = depth;
        return;
    }
    const size_t split = std::bit_floor(drops.size() - 1);
    FillDrops(drops.first(split), depth - 1);
    FillDrops(drops.subspan(split), depth - 1);
}

Row Merge(const Setup& s, const Row& low, const Row& high, bool lazy) {
    Row result(low.size() * high.size());
#pragma omp parallel for
    for (size_t i = 0; i < result.size(); ++i) {
        auto product = s.cc->EvalMultNoRelinNoCheck(low[i % low.size()], high[i / low.size()]);
        if (!lazy)
            s.cc->RelinearizeInPlace(product);
        s.cc->ModReduceInPlace(product);
        result[i] = std::move(product);
    }
    return result;
}

Row Build(const Setup& s, std::span<Row> rows, bool lazy) {
    if (rows.size() == 1)
        return std::move(rows[0]);
    const size_t split = std::bit_floor(rows.size() - 1);
    auto low = Build(s, rows.first(split), false);
    auto high = Build(s, rows.subspan(split), false);
    return Merge(s, low, high, lazy);
}

Ciphertext WeightedSum(const Setup& s, const Row& row, std::span<const int64_t> lut, size_t offset,
                       size_t stride) {
    auto result = s.cc->GetScheme()->MultByInteger(row[0], lut[offset]);
    for (size_t i = 1; i < row.size(); ++i)
        s.cc->EvalAddInPlaceNoCheck(result,
                                    s.cc->GetScheme()->MultByInteger(row[i], lut[offset + i * stride]));
    return result;
}
void ValidateTree(const Setup& s, const SelectorMatrix& selectors, uint64_t entryCount) {
    TreeLevelDrops(selectors.size());
    uint64_t entries = 1;
    Ciphertext first;
    for (const auto& row : selectors) {
        if (row.empty() || row.size() > (uint64_t{1} << 32) / entries)
            throw std::invalid_argument("Empty selector row or tree larger than 32 bits");
        entries *= row.size();
        for (const auto& ct : row) {
            if (!ct || ct->GetCryptoContext() != s.cc || ct->GetKeyTag() != s.publicKey->GetKeyTag() ||
                ct->GetElements().size() != 2 || ct->GetNoiseScaleDeg() != 1)
                throw std::invalid_argument("Invalid tree selector");
            if (!first)
                first = ct;
            if (ct->GetLevel() != first->GetLevel() || ct->GetScalingFactor() != first->GetScalingFactor())
                throw std::invalid_argument("Tree selectors must have a common level and scale");
        }
    }
    if (entries != entryCount)
        throw std::invalid_argument("LUT size does not match the selector rows");
}

void DropTreeLevels(const Setup& s, SelectorMatrix& selectors) {
    const auto drops = TreeLevelDrops(selectors.size());
    for (size_t digit = 0; digit < selectors.size(); ++digit) {
        if (drops[digit] == 0)
            continue;
        for (auto& ct : selectors[digit])
            ct = s.cc->LevelReduce(ct, nullptr, drops[digit]);
    }
}

template <class Entry>
Ciphertext EncodedSum(const Setup& s, const Row& row, std::span<const Entry> lut, size_t offset,
                      size_t stride) {
    Ciphertext result;
    for (size_t i = 0; i < row.size(); ++i) {
        Ciphertext product;
        if constexpr (std::is_same_v<Entry, lbcrypto::Plaintext>)
            product = s.cc->EvalMult(row[i], lut[offset + i * stride]);
        else
            product = s.cc->EvalMultNoRelinNoCheck(row[i], lut[offset + i * stride]);
        if (result)
            s.cc->EvalAddInPlaceNoCheck(result, product);
        else
            result = std::move(product);
    }
    if constexpr (std::is_same_v<Entry, Ciphertext>)
        s.cc->RelinearizeInPlace(result);
    s.cc->ModReduceInPlace(result);
    return result;
}

template <class Entry>
Ciphertext EncodedTree(const Setup& s, SelectorMatrix selectors, std::span<const Entry> lut) {
    ValidateTree(s, selectors, lut.size());
    const auto branchDepth = selectors.size() == 1 ? 0 : std::bit_width(selectors.size() - 1) - 1;
    if (selectors[0][0]->GetLevel() + branchDepth != detail::LutEntryLevel(s))
        throw std::invalid_argument("Selector level does not match the encoded LUT entries");
    DropTreeLevels(s, selectors);
    if (selectors.size() == 1)
        return EncodedSum(s, selectors[0], lut, 0, 1);
    auto rows = std::span(selectors);
    const size_t split = (selectors.size() + 1) / 2;
    auto low = Build(s, rows.first(split), true);
    auto high = Build(s, rows.subspan(split), false);
    // Entry multiplication consumes a level on the high branch.
    for (auto& ct : low)
        ct = s.cc->LevelReduce(ct, nullptr, 1);
    Ciphertext result;
#pragma omp parallel
    {
        Ciphertext partial;
#pragma omp for schedule(static)
        for (size_t i = 0; i < low.size(); ++i) {
            auto weighted = EncodedSum(s, high, lut, i, low.size());
            auto product = s.cc->EvalMultNoRelinNoCheck(low[i], weighted);
            if (partial)
                s.cc->EvalAddInPlaceNoCheck(partial, product);
            else
                partial = std::move(product);
        }
#pragma omp critical
        {
            if (partial) {
                if (result)
                    s.cc->EvalAddInPlaceNoCheck(result, partial);
                else
                    result = std::move(partial);
            }
        }
    }
    s.cc->RelinearizeInPlace(result);
    s.cc->ModReduceInPlace(result);
    return result;
}
} // namespace

std::vector<uint32_t> TreeLevelDrops(uint32_t digits) {
    if (digits < 1 || digits > 32)
        throw std::invalid_argument("A tree requires 1..32 digits");
    std::vector<uint32_t> drops(digits);
    if (digits > 1) {
        const size_t split = (digits + 1) / 2;
        const uint32_t sideDepth = std::bit_width(digits - 1) - 1;
        FillDrops(std::span(drops).first(split), sideDepth);
        FillDrops(std::span(drops).subspan(split), sideDepth);
    }
    return drops;
}

Ciphertext EvaluateCleartextTree(const Setup& s, SelectorMatrix selectors, std::span<const int64_t> lut) {
    ValidateTree(s, selectors, lut.size());
    for (auto value : lut)
        if (value < 0)
            throw std::invalid_argument("Cleartext tree entries must be nonnegative");

    DropTreeLevels(s, selectors);
    if (selectors.size() == 1)
        return WeightedSum(s, selectors[0], lut, 0, 1);

    const size_t split = (selectors.size() + 1) / 2;
    auto rows = std::span(selectors);
    // Keep the low side unrelinearized until the final contraction. The high
    // side has two components, so the final products have degree at most three.
    auto low = Build(s, rows.first(split), true);
    auto high = Build(s, rows.subspan(split), false);
    Ciphertext result;
#pragma omp parallel
    {
        Ciphertext partial;
#pragma omp for schedule(static)
        for (size_t i = 0; i < low.size(); ++i) {
            auto weighted = WeightedSum(s, high, lut, i, low.size());
            auto product = s.cc->EvalMultNoRelinNoCheck(low[i], weighted);
            if (partial)
                s.cc->EvalAddInPlaceNoCheck(partial, product);
            else
                partial = std::move(product);
        }
#pragma omp critical
        {
            if (partial) {
                if (result)
                    s.cc->EvalAddInPlaceNoCheck(result, partial);
                else
                    result = std::move(partial);
            }
        }
    }
    // Contract with the LUT before relinearizing, including the two-digit case.
    s.cc->RelinearizeInPlace(result);
    s.cc->ModReduceInPlace(result);
    return result;
}
namespace detail {
uint32_t LutEntryLevel(const Setup& s) {
    return s.depth - s.parameters.levelBudget[1] - (s.digitWidths.size() == 1 ? 1 : 2);
}

Ciphertext EvaluatePlaintextTree(const Setup& s, SelectorMatrix selectors,
                                 std::span<const lbcrypto::Plaintext> lut) {
    return EncodedTree(s, std::move(selectors), lut);
}

Ciphertext EvaluateCiphertextTree(const Setup& s, SelectorMatrix selectors, std::span<const Ciphertext> lut) {
    return EncodedTree(s, std::move(selectors), lut);
}
} // namespace detail
} // namespace large_lut
