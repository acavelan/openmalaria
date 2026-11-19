# Rewrite KPI ledger (plan v3.3)

One row per phase (contract: `openmalaria-rewrite` skill). Three KPIs per row.

- **model LOC** = `cloc --quiet model --include-lang="C++,C/C++ Header" | grep '^SUM'` col 5 — code lines only (comments/blanks never count; `// parity:` annotations are LOC-free). Hard cap Σ ≤ 14.0k incl. `model/xml/` (xml/ ≤ 2.6k). `unittest/` LOC tracked beside it — model logic relocated into tests counts as a breach.
- **CPU-s** = user CPU seconds (`/usr/bin/time -f %U`) of the single benchmark: `util/example/example_scenario.xml` (popSize 2000, iseed 1 — built-in values, scenario not to be edited; covers dynamic EIR, 3 species, ExactNv0, ITN deployments = full hotspot map), `taskset`-pinned, governor/turbo off, median-of-5, warmup/main split, **build A2** (PGO; A0 = plain `-O3 -DNDEBUG -ffp-contract=off` is the T1 gate build). Never ctest wall. Noise floor = max(1.0%, 2×MAD); below it = revert. The benchmark scenario also runs bit-exact-diff in the parity runner (timing is void if its output diverges).
- **Readability** = qualitative note (new free-function equations, `// parity:` anchors kept, no statics).

Build-recipe hash (compiler/flags/lib versions/source commit of the T1 reference binary) is stamped on every row; hash mismatch invalidates parity claims.

| Phase | Description | model LOC | unittest LOC | CPU-s (A2, example_scenario w/m) | vs A2 prev | Readability note | Verdict |
|---|---|---|---|---|---|---|---|
| 0 | baseline, build A0 `-O3 -DNDEBUG -ffp-contract=off` (recipe `fb6a8d73c3c0bd90`, ref bin `/tmp/opencode/om-ref/A0/`, commit cfc3c8e0 clean) | 18,813 | 3,529 | A0 12.33 (MAD .010) | 1.000 | reference | reference |
| 0a | flag lever A1 = A0+LTO — T1 bit-check **57/57 pass** | — | — | 12.03 (MAD .010) | **A1/A0 = 1.025** | — | kept (free) |
| 0b | flag lever A2 = A1+PGO (corpus: ctest 51 + unittest + example hot-run) — T1 bit-check **57/57 pass** ⇒ **A2 is the timing yardstick** (recipe `52a4b6f5b1ac0649`) | — | — | 11.34 (MAD .010) | **A2/A0 = 1.087** | — | kept; note: slightly under the 1.10 hope |
| 0c | branch-miss pricing (`perf stat`, example@2000): A0 120.8M misses vs A2 120.3M flat; IPC 3.56→3.63; PGO win = −5.8 % instructions, not miss-rate. Residual miss ceiling ≈ 0.33 s/12 s (**≤3 %**) — de-virt (phase 4) must win via instruction/indirect-call elimination + IPO, *not* branch prediction; KPI 1.15–1.35 pre-break already reflects this | — | — | — | — | — | de-virt priced |
| 0d | harness: `util/parity.py` (T1 two-tier runner, recipe-hash guard), `util/toleranceTable.py` → `test/tolerance_example_p2000.json` (189,800 keys, CV med 2.7 %/p95 316 % ⇒ blanket 2e-2 falsified), `--seed` CLI (bit-proven ≡ XML iseed), EIR-uninit UB fix (**value-neutral, zero re-bless**), single-run gate mode: **53/53 bin-identical**; ctest noise 9.36 s | 18,839 (+26 seed/UB-fix) | grows with T2 work | — | — | +`--seed`, −UB | green |

Known: `12a` continuous-vs-resumed golden divergence (22 measures, ≤8 %; checkpoint veto covers it; re-bless scheduled at phase 8).
