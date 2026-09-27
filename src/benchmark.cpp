// SPDX-License-Identifier: BSD-2-Clause
#include "lut.hpp"
#include <algorithm>
#include <charconv>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <omp.h>
#include <random>
#include <stdexcept>
#include <string>
#include <string_view>

using namespace large_lut;

namespace {
using Clock = std::chrono::steady_clock;

struct Options {
    Parameters parameters;
    uint32_t runs = 1;
    uint64_t seed = 1;
    int threads = 0;
    bool slotVarying = false;
    bool help = false;
};

void Help() {
    std::cout << R"(Usage: lut-bench [options]

Generate a random LUT, evaluate RLWE inputs, and verify decrypted outputs.
Defaults are toy parameters with no security guarantee.

  --input-bits B       LUT has 2^B entries, B=1..32 (default: 8)
  --digit-bits B       Digit width, B=1..8 (default: 4; partial top digit supported)
  --kind KIND          cleartext, plaintext, or ciphertext (default: cleartext)
  --secure             Check security bounds; defaults to log-n=16, large-digits=4
  --log-n L            Ring dimension 2^L, L=8..30 (toy default: 8)
  --runs R             Number of verified evaluations (default: 1)
  --seed S             Random LUT and input seed (default: 1)
  --threads T          OpenMP threads (default: OpenMP environment/runtime)
  --h H                Secret-key Hamming weight (default: 192)
  --cleaning C         Cleaning iterations, C=0..8 (-1: automatic, default)
  --scale-bits B       Scaling modulus bits, B=40..59 (0: automatic, default)
  --large-digits D     Key-switch decomposition digits (toy default: 3)
  --level-budget E D   Encoding and decoding level budgets (default: automatic)
  --slot-varying       Offset encoded LUT values by coefficient index
  --help              Show this help

Each run verifies N inputs. Domains up to N are covered completely; larger
domains use random inputs and carry boundaries. Evaluation timing excludes
setup, LUT preparation, input encryption, and output decryption/verification.
)";
}

template <class T> T Number(std::string_view text, std::string_view option) {
    T value{};
    auto [end, error] = std::from_chars(text.data(), text.data() + text.size(), value);
    if (error != std::errc{} || end != text.data() + text.size())
        throw std::invalid_argument("Invalid value for " + std::string(option) + ": " + std::string(text));
    return value;
}

Options Parse(int argc, char** argv) {
    Options o;
    bool ringGiven = false, largeDigitsGiven = false;
    for (int i = 1; i < argc; ++i) {
        const std::string_view option = argv[i];
        auto next = [&]() -> std::string_view {
            if (++i >= argc)
                throw std::invalid_argument("Missing value for " + std::string(option));
            return argv[i];
        };
        if (option == "--help") {
            o.help = true;
            return o;
        } else if (option == "--input-bits")
            o.parameters.inputBits = Number<uint32_t>(next(), option);
        else if (option == "--digit-bits")
            o.parameters.digitBits = Number<uint32_t>(next(), option);
        else if (option == "--kind") {
            const auto kind = next();
            if (kind == "cleartext")
                o.parameters.lutKind = LutKind::Cleartext;
            else if (kind == "plaintext")
                o.parameters.lutKind = LutKind::Plaintext;
            else if (kind == "ciphertext")
                o.parameters.lutKind = LutKind::Ciphertext;
            else
                throw std::invalid_argument("Expected cleartext, plaintext, or ciphertext for --kind");
        } else if (option == "--secure")
            o.parameters.toy = false;
        else if (option == "--log-n") {
            const auto logN = Number<uint32_t>(next(), option);
            if (logN < 8 || logN > 30)
                throw std::invalid_argument("--log-n must be in 8..30");
            o.parameters.ringDim = uint32_t{1} << logN;
            ringGiven = true;
        } else if (option == "--runs") {
            o.runs = Number<uint32_t>(next(), option);
            if (o.runs == 0)
                throw std::invalid_argument("--runs must be positive");
        } else if (option == "--seed")
            o.seed = Number<uint64_t>(next(), option);
        else if (option == "--threads") {
            o.threads = Number<int>(next(), option);
            if (o.threads < 1)
                throw std::invalid_argument("--threads must be positive");
        } else if (option == "--h")
            o.parameters.h = Number<uint32_t>(next(), option);
        else if (option == "--cleaning")
            o.parameters.cleaning = Number<int>(next(), option);
        else if (option == "--scale-bits")
            o.parameters.scaleBits = Number<uint32_t>(next(), option);
        else if (option == "--large-digits") {
            o.parameters.largeDigits = Number<uint32_t>(next(), option);
            largeDigitsGiven = true;
        } else if (option == "--level-budget") {
            const auto encoding = Number<uint32_t>(next(), option);
            const auto decoding = Number<uint32_t>(next(), option);
            o.parameters.levelBudget = {encoding, decoding};
        } else if (option == "--slot-varying")
            o.slotVarying = true;
        else
            throw std::invalid_argument("Unknown option: " + std::string(option));
    }
    if (!o.parameters.toy) {
        if (!ringGiven)
            o.parameters.ringDim = 1u << 16;
        if (!largeDigitsGiven)
            o.parameters.largeDigits = 4;
    }
    if (o.slotVarying && o.parameters.lutKind == LutKind::Cleartext)
        throw std::invalid_argument("--slot-varying requires a plaintext or ciphertext LUT");
    ValidateParameters(o.parameters);
    return o;
}

double Seconds(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}

const char* KindName(LutKind kind) {
    switch (kind) {
    case LutKind::Cleartext:
        return "cleartext";
    case LutKind::Plaintext:
        return "plaintext";
    case LutKind::Ciphertext:
        return "ciphertext";
    }
    throw std::invalid_argument("Invalid LUT kind");
}

template <class Lut>
void Benchmark(const Setup& s, const Options& o, const std::vector<int64_t>& table, const Lut& lut,
               double preparationSeconds) {
    std::cout << "lut_preparation_s=" << preparationSeconds << '\n' << std::flush;
    auto evaluator = s;
    evaluator.secretKey.reset();
    std::mt19937_64 random(o.seed);
    std::vector<double> times;
    for (uint32_t run = 0; run < o.runs; ++run) {
        std::vector<int64_t> values(s.parameters.ringDim);
        for (size_t i = 0; i < values.size(); ++i)
            values[i] = s.lutSize <= values.size() ? i % s.lutSize : random() % s.lutSize;
        if (s.lutSize > values.size()) {
            size_t i = 0;
            values[i++] = 0;
            values[i++] = s.lutSize - 1;
            for (uint32_t bit = 0; bit < s.parameters.inputBits; ++bit) {
                values[i++] = (int64_t{1} << bit) - 1;
                values[i++] = int64_t{1} << bit;
            }
        }
        std::shuffle(values.begin(), values.end(), random);
        auto input = EncryptInput(s, values);
        const auto start = Clock::now();
        auto output = Evaluate(evaluator, input, lut);
        const double elapsed = Seconds(start);
        auto actual = DecryptOutput(s, output);
        if (actual.size() != values.size())
            throw std::runtime_error("Unexpected output coefficient count");
        for (size_t i = 0; i < values.size(); ++i) {
            const auto expected = (table[values[i]] + (o.slotVarying ? i : 0)) % s.lutSize;
            if (uint64_t(actual[i]) != expected)
                throw std::runtime_error("Run " + std::to_string(uint64_t(run) + 1) + ", coefficient " +
                                         std::to_string(i) + ": expected " + std::to_string(expected) +
                                         ", got " + std::to_string(actual[i]));
        }
        times.push_back(elapsed);
        std::cout << "run=" << uint64_t(run) + 1 << " evaluation_s=" << elapsed
                  << " checked=" << values.size() << " mismatches=0\n"
                  << std::flush;
    }
    std::sort(times.begin(), times.end());
    const double median = (times[(times.size() - 1) / 2] + times[times.size() / 2]) / 2;
    std::cout << "median_s=" << median
              << " mean_s=" << std::accumulate(times.begin(), times.end(), 0.0) / times.size()
              << " min_s=" << times.front() << " max_s=" << times.back()
              << " lookups_per_second=" << s.parameters.ringDim / median
              << " checked=" << uint64_t(o.runs) * s.parameters.ringDim << " mismatches=0\n";
}
} // namespace

int main(int argc, char** argv) {
    try {
        const auto o = Parse(argc, argv);
        if (o.help) {
            Help();
            return 0;
        }
        if (o.threads) {
            omp_set_dynamic(0);
            omp_set_num_threads(o.threads);
        }
        std::cout << std::fixed << std::setprecision(6)
                  << "mode=" << (o.parameters.toy ? "toy (insecure)" : "secure")
                  << " kind=" << KindName(o.parameters.lutKind) << " input_bits=" << o.parameters.inputBits
                  << " digit_bits=" << o.parameters.digitBits << " ring_dim=" << o.parameters.ringDim
                  << " threads=" << omp_get_max_threads() << " runs=" << o.runs << " seed=" << o.seed
                  << " slot_varying=" << o.slotVarying << '\n'
                  << std::flush;
        auto start = Clock::now();
        const auto s = CreateSetup(o.parameters);
        std::cout << "setup_s=" << Seconds(start) << " entries=" << s.lutSize
                  << " scale_bits=" << s.parameters.scaleBits << " cleaning=" << s.cleaning
                  << " h=" << s.parameters.h << " large_digits=" << s.parameters.largeDigits
                  << " level_budget=" << s.parameters.levelBudget[0] << ',' << s.parameters.levelBudget[1]
                  << " depth=" << s.depth << " logQP=" << s.logQP << '\n'
                  << std::flush;
        start = Clock::now();
        std::mt19937_64 random(o.seed);
        std::vector<int64_t> table(s.lutSize);
        for (auto& value : table)
            value = random() % s.lutSize;
        if (o.parameters.lutKind == LutKind::Cleartext) {
            Benchmark(s, o, table, table, Seconds(start));
        } else {
            const auto encode = [&] {
                return EncodeLut(s, [&](uint64_t index, uint32_t coefficient) -> int64_t {
                    return (table[index] + (o.slotVarying ? coefficient : 0)) % s.lutSize;
                });
            };
            if (o.parameters.lutKind == LutKind::Plaintext) {
                const auto lut = encode();
                Benchmark(s, o, table, lut, Seconds(start));
            } else {
                const auto lut = EncryptLut(s, encode());
                Benchmark(s, o, table, lut, Seconds(start));
            }
        }
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}
