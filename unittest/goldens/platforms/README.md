# Bit-exact platform baselines

**Zero tolerance.** Each build selects one frozen reference before running tests.
Every RNG value and state word, inverse-CDF spot, availability percentile and CFR
fixture is compared bit-for-bit. The extra live model-versus-adapter checks remain.

| CMake selection | Architecture | GSL | Compiler at capture | RNG values differing from original | Percentiles differing |
|---|---|---|---|---:|---:|
| `ubuntu-22.04` | x86-64 | 2.7.1 | GCC 11.4 | 2 | 0 |
| `ubuntu-24.04` | x86-64 | 2.7.1 | GCC 13.3 | 0 | 0 |
| `ubuntu-24.04-arm` | AArch64 | 2.7.1 | GCC 13.3 | 1,051 | 5 |
| `macos-15` | ARM64 | 2.8 | AppleClang 17 | 1,122 | 5 |
| `macos-26` | ARM64 | 2.8 | AppleClang 21 | 1,122 | 5 |

Use `cmake ... -DOM_GOLDEN_PLATFORM=ubuntu-22.04`, for example. CI sets
`OM_GOLDEN_PLATFORM` to its matrix runner label. An empty selection retains the
original Linux x86-64/GSL 2.8 fixtures. Unknown names fail configuration; a library
update that changes results fails the exact tests. There is no result-based
selection, tolerance, runtime reference generation, or automatic re-blessing.

## Provenance

[Native capture run](https://github.com/acavelan/openmalaria/actions/runs/45777510549),
driver commit `057172195d8f06b55b8e7968be7487ecdf7d9895`:

- Original model/engines: `cfc3c8e0caee6fdcd6b488a5ab50fdbf2b7a5a90`.
- Original RNG adapter/generator: `6e24e51eb576a12fd701ab54a1b6c5dd375ca081`,
  extracted from git, including its unchanged distribution schedule.
- Percentiles: `percentiles_gen.cpp` calls the original model's
  `PerHostAnophParams::calcAvailabilityPercentiles()` directly.
- Flags: `-O3 -DNDEBUG -ffp-contract=off`, no LTO/PGO.
- Each `provenance.txt` records runner image, OS/CPU, compiler, linked libraries,
  source hashes, complete capture hashes and stored fixture hashes. The original
  RNG generator's CSV header has historical x86/GCC wording; the provenance file
  records the actual platform/compiler.

All 14,012 rows have the **same schedule and state words** as the original CSV.
`rng_values.csv` is a sparse, lossless encoding of the native capture: only changed
values are stored, keyed by distribution, parameters and draw index. Configuration
expands it using `apply-platform.py`; the full original state table is retained.
Decimal strings have 17 significant digits and round-trip to the captured doubles.
The reconstructed data rows were verified against every complete native capture.
Uniform/discrete samples and CFR use the original expectations on every platform.

`availability_percentiles.inc` contains the complete 100-entry native table.
The original 100 literals were moved, unchanged, to `../availability_percentiles.inc`.
The pinned-source generator also reproduces all original values on the original
local toolchain.

## Reproduce an independent capture

With dependencies installed and both pinned commits available locally:

```sh
bash unittest/golden-gen/build-goldens.sh /tmp/opencode/new-baseline
```

Use a new output directory. The script builds extracted original sources, records
provenance, and produces both complete captures and the compact fixture files.
It never overwrites the checked-in references. Updating a baseline requires an
explicitly reviewed fixture change; building/testing only reads existing fixtures.
