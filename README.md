Evaluating Larger Lookup Tables using CKKS
---

Proof of concept implementation of the paper [Evaluating Larger Lookup Tables using CKKS](https://tches.iacr.org/index.php/TCHES/article/view/12693).
It is based on the more recent [OpenFHE v1.6.0](https://github.com/openfheorg/openfhe-development/releases/tag/v1.6.0), and hence has better performance than in the paper.
Integration in the main OpenFHE library is planned.

## Installation

With OpenFHE v1.6.0 installed systemwide (64-bit native integers), Clang with C++20 support, CMake, and OpenMP:

```sh
cmake -S . -B build -DCMAKE_CXX_COMPILER=clang++
cmake --build build -j
```

Release and native CPU optimizations are enabled by default.

## Running benchmarks

```sh
# Quick toy run: N=256, no security guarantee
./build/lut-bench --input-bits 8 --digit-bits 4 --threads 16

# Secure 16-bit random LUT: N=65536, three verified runs
./build/lut-bench --secure --input-bits 16 --digit-bits 4 --threads 16 --runs 3

# Encoded LUTs (also supports --kind ciphertext)
./build/lut-bench --kind plaintext --input-bits 5 --digit-bits 2 --slot-varying --threads 16

./build/lut-bench --help
```

Supports 1–32 input bits, 1–8 bits per digit, and cleartext, plaintext, or ciphertext LUTs.
Toy mode is the default.
Use `--secure` explicitly for secure runs.
Large LUTs can require substantial RAM, especially with encoded entries.

Each run verifies the decrypted outputs.
Domains no larger than N are checked in full.
Larger domains use random inputs and carry boundaries.
`--seed` selects the random LUT and inputs.
Timings report setup and LUT preparation separately from RLWE-to-RLWE evaluation, which excludes input encryption and output verification.
Without `--threads`, OpenMP uses its environment/runtime settings.

Three small toy checks are available with `OMP_NUM_THREADS=16 ctest --test-dir build --output-on-failure`.

## Docker

The image builds OpenFHE v1.6.0 and the benchmark executable:

```sh
docker build -t large-lut .
docker run --rm large-lut --input-bits 8 --digit-bits 4 --threads 16
docker run --rm large-lut --secure --input-bits 16 --digit-bits 4 --threads 16 --runs 3
```

The container accepts the same CLI options as the local executable.

## Cite

```bibtex
@article{LargeLutCKKS,
	title        = {Evaluating Larger Lookup Tables using CKKS},
	author       = {Dumezy, Jules and Alexandru, Andreea and Polyakov, Yuriy and Clet, Pierre-Emmanuel and Chakraborty, Olive and Boudguiga, Aymen},
	year         = 2026,
	month        = {Jan.},
	journal      = {IACR Transactions on Cryptographic Hardware and Embedded Systems},
	volume       = 2026,
	number       = 1,
	pages        = {559–591},
	doi          = {10.46586/tches.v2026.i1.559-591},
	url          = {https://tches.iacr.org/index.php/TCHES/article/view/12693},
}
```
