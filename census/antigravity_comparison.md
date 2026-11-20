# antigravity-refactor (Rust) vs. this rewrite — comparison memo (phase 7c-A)

Source audited read-only: `/home/razlock/git/blake-antigravity` @ `antigravity-refactor`
(telemetry ledger `migration_telemetry.json`; 12 `#![forbid(unsafe_code)]` crates + PyO3).
Their claim: **1.46× faster, 67.7 % leaner RSS, 201/201 parity, 0.00e+00 error**.
Question for us: which of their levers are adoptable into the C++ rewrite, and what is
each worth **on like-for-like, bit-exact terms**?

## 1. The two staircases — NOT directly comparable

Their numbers (telemetry `benchmarks[]`: **wall-time** of a **5-scenario** mix —
VecFullTest, 2ITNs, Molineaux, Vivax, EffectiveDrug — median of 3, whole process):

| Stage | What | wall ms | RSS MB | speedup |
|---|---|---|---|---|
| S1 | legacy C++ `-O3` (CMake default = **generic x86-64 SSE2** baseline) | 6241 | 17.41 | 1.00 |
| S2 | Safe Rust AoS, multi-CGU, no LTO | 5917 | 6.33 | 1.05 |
| S3 | Rust SoA gather-scratch + fat-LTO + CGU=1 + panic=abort | 5188 | 5.78 | 1.20 |
| S4 | S3 + PGO **+ `-Ctarget-cpu=native` (AVX2/FMA)** | 4262 | 5.62 | **1.46** |

Ours (PERF.md: **user CPU-s**, single KPI benchmark `example_scenario.xml` @ popSize 2000,
`taskset -c 5`, median-of-5; ratio vs frozen baseline-A0 12.33):

| Stage | What | CPU-s | ratio vs baseline-A0 |
|---|---|---|---|
| A0 | baseline tree `-O3 -DNDEBUG -ffp-contract=off` | 12.33 | 1.000 |
| A1 | baseline + `-flto=auto` | 12.03 | 1.025 |
| A2 | baseline + PGO (frozen yardstick) | 11.34 | 1.087 |
| p7b-green A0-class | current rewritten tree, gate flags | 12.47 | 0.989 |
| p7b-green "A2" | current tree, a2-recipe | **11.98** | 1.029 (0.947 vs frozen A2 11.34) |

Why the bars don't line up:

1. **Metric**: they sum wall-time over 5 *small* scenarios — process startup + scenario
   XML parse is a large fixed share there. The Rust side replaces Xerces+xsdcxx with a
   light pure-Rust parser, so S2..S4 inherit a startup/parse win that has nothing to do
   with the simulation loop. We measure **user CPU-s of one pop-2000 run**, startup/parse
   amortized to ~zero. A wall-mix ratio and a hot-loop CPU-s ratio are different quantities.
2. **ISA asymmetry on their baseline** (their S1 vs S4): C++ built for generic x86-64
   (SSE2), Rust built `-Ctarget-cpu=native` (AVX2+FMA+tuning). A pure compiler/ISA effect
   is charged to the port. We price this exactly, same-tree ±flags, in 7c-B (§4 below).
3. **Their S3→S4 step bundles PGO with native-ISA** (1.20→1.46 = 1.22× combined) — the
   telemetry does not separate the two ingredients. Our experiment separates them.
4. **Parity bar**: they compare **survey output text within 1e-6** (201/201 "0.00e+00");
   we require **bit-identical `output.bin`**. Their engine is *not* the same computation
   (§3 flags): hand-LU solver, a reordered product, hand-rolled CDFs, value-"fixed" latent
   bugs. Their measured speed is of a *different program*; ours must survive `cmp`.
5. **Baseline tree ≠ current tree**: their S1 is the same commit our 12.33 is, but our tree
   has since deleted the polymorphism/statics/if-ladders both fat-LTO/PGO (their side) and
   devirtualization (ours) would otherwise mine. Their levers were priced against a
   baseline that still had that goldmine; ours no longer does.

## 2. Per-lever adoption table

| Lever (their stage) | Their evidence | Our status | Verdict | Bit-exact? |
|---|---|---|---|---|
| Two-pass SoA scratch sweeps: flat contiguous `[host]` and `[host×genotype]` f64 buffers staged then reduced; zero hot allocs (S2→S3, ~half of 1.143) | 1.05→1.20 incl. fat-LTO | B1a did the alloc half (WH `y_lag` scratch reuse, −2.4 % KEEP); `updateKappa`/`vectorUpdate` already scratch-shaped; **7c-C completes the two-pass loop shape + makes it vectorize** | **ADOPT (loop shape only)** | yes — elementwise staging + scalar ordered reduce, no reassociation |
| Fat-LTO + CGU=1 + panic=abort (S3) | in 1.05→1.20 | GCC `-flto=auto` already in A1/A2 (worth 1.025 on baseline tree, ~1.008 measured on current tree); Rust CGU=1 is LLVM's way of *enabling* cross-TU work — no separate GCC equivalent to buy | **DONE** | yes (LTO passes T1 on both trees) |
| PGO (S4, bundled w/ native) | in 1.20→1.46 | **Audit finding (7c, this phase): our p7b "A2" yardstick binary never loaded the profile data** — `-fprofile-use` looked for `.gcda` beside `build-p7b-a2` objects, none existed (profiles sat path-encoded for `build-p7b-pgogen` in `/tmp/opencode/om-pgo`), `-Wno-missing-profile` silenced it. "PGO collapsed to 1.003" was **never actually tested on the current tree**. 7c re-prices PGO with correctly wired profiles (rename + `-fprofile-dir`, mechanism probe-verified) | **REPRICE (7c-E)** | must pass T1 bit-check per build — PGO passed on baseline, re-prove on current |
| `-Ctarget-cpu=native` ISA edge (S4) | in 1.20→1.46, inseparable | 7c-B builds current tree a2-recipe ± `-march=native -mtune=native -ffp-contract=off`, full T1 a+b, bench. Caveats priced: our `dgemm` (~20 %) runs in **runtime-dispatching BLAS (already AVX2)** and `exp/pow` in **glibc IFUNC libm (already AVX2)** — native mainly affects model/ loops, so our ISA delta should be SMALLER than the headline suggests | **MEASURE, propose USER-class "A2n" only with user ratification (plan §2.7 bans `-march` from GATE binaries; contract stays `-ffp-contract=off` + FMA-emission zero check)** | conditional — only if T1 63/63 on the native build |
| Hand-LU `Nv0` solve replacing dgemm (their solver rewrite) | part of parity-debugged port | our solver flip is the **phase-8 sanctioned break** (F3 `solve_Nv0_linear`, ≈21 % of CPU → ≈1.19 alone) | **DEFER to phase 8** | no — different computation |
| Hand-rolled CDF grids (`ugaussian_P`, lognormal/gamma/beta ports) | ride the 1e-6 bar | phase 9/T3: inline ports allowed **only where bit-parity grids pass** (plan F2), else keep GSL | **DEFER to phase 9** | no |
| Sampler replacements (KS-level) | ride the 1e-6 bar | phase 8/9 territory (draw order frozen pre-break) | **DEFER** | no |
| tbv multiply reorder | invisible at 1e-6 | one FP op-order change = ours diverges chaotically | **REJECT pre-break** | no |
| Structural statics→context, de-polymorph, flat interventions | "1188 hazards", their crate split | we did it our way (phases 3–7, RunContext + enum-Kind + flat tables) — zero polymorphism left in model/, verified independently of them | **DONE (our route)** | yes |
| RSS −68 % (5.6 MB) | telemetry | byproduct of small-scenario mix + Rust allocator + their scratch reuse; not a KPI for us (RSS ≈ baseline parity since p5); slot-traversal (phase 8) removes the residual churn | **NOT ADOPTED (no KPI)** | — |

## 3. The four "only works at 1e-6" flags

Each of these is invisible under their text-parity bar but **breaks bit-exactness**
against the reference binary; none may enter our tree before the sanctioned break:

1. **Hand-LU Nv0 solver** — replaces the dgemm/gesv path with a different elimination
   order → different roundings on the exactness-critical linear solve (their 380-magnitude
   iteration-1 divergences were exactly this class). Phase 8 (F3) owns our flip.
2. **tbv multiply reorder** — one product computed in a different order; ~1 ULP per
   application, amplified by the model's feedback loops. Phase 8 reassociation.
3. **Hand-rolled CDFs** (inverse-normal etc.) — a different polynomial than GSL → ULP
   deltas on every draw they gate. Only if phase-9 T3 grids show bit-parity; else GSL stays.
4. **Value-changing latent-bug fixes** (8 per their ch1: uninitialized EIR accumulators,
   skipped post-warmup reset, stale cohort masks, resample spin) — under 1e-6 text parity a
   "fix" is invisible only if the bug's numeric effect is below 1e-6; several of theirs are
   not. Ours: every latent-bug fix is value-neutral or deferred to the phase-8 re-bless
   (our EIR-uninit fix was proven value-neutral, zero re-bless).

## 4. Like-for-like pricing of the ISA edge (7c-B result)

Same tree, same a2-recipe flags, only `± -march=native -mtune=native` (with
`-ffp-contract=off` — `objdump -d` of the whole binary: **0** FMA/FMADD/FMSUB emissions):

| Build | median CPU-s (of 5, MAD) | T1 a+b | Δ |
|---|---|---|---|
| `build-p7c-a2` (no -march) | **11.98** (.010; raw 11.96–12.00) | 63/63 (yardstick proven at p7b) | baseline |
| `build-p7c-a2n` (`-march=native`) | **12.12** (.010; raw 12.10–12.13) | **63/63 bit-identical** | **+1.17 % SLOWER** |

The distributions are disjoint (a2 max 12.00 < a2n min 12.10). `perf stat`: a2n retires
**1.05 % fewer instructions** (221.7 B vs 224.0 B) yet takes **1.28 % more cycles**
(IPC 3.53 → 3.45) — the AVX2/alderlake codegen shrank the code but scheduled it worse for
this scalar-FP dependency-bound model.

**Verdict: not proposed for adoption.** The plan's `-march`/`-mtune` ban costs nothing on
this tree; there is no "A2n" to ratify. This directly prices the biggest single
ingredient of their S3→S4 (1.22×) step: the same ISA class applied to *our* current tree
is a penalty, not a gain. Consistent with §2's caveat — the two hottest kernels
(runtime-dispatching AVX2 BLAS ≈20 %, glibc IFUNC libm `exp/pow`) are already native-speed
on **both** sides of the A/B, and IPC 3.53 leaves vectorization no room: the remaining
time is scalar chains + memory, where AVX2 mainly adds store-forwarding friction.
What remains of their bundled step after removing the ISA edge is PGO (re-priced honestly
in 7c-E after fixing our profile wiring) and whatever the small-scenario wall-mix effect
contributes there.

