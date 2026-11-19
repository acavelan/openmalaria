# Phase 5 HANDOFF — cycle breaks + Infection perf fix (branch `p5-cycles` off `p4b-green`)

Crash-recovery file. Update EVERY green commit. Baseline p4b-green: model **18,599** / unittest 3,760; A0 bench 12.33s; p4b bench 12.80s (session re-measure `p5-before` = **12.77s** mad 0.000). Full phase-4b handoff (per-hierarchy lessons): `git show p4b-green:HANDOFF.md`.

## Environment
- Build: `build-p5/` = Release `-ffp-contract=off` `-DOM_BOXTEST_BINARY=ON` (configured).
- unittest MUST run from `build-p5/unittest/` (csv resources there) → `./unittest` (91 tests).
- ctest: `cd build-p5 && ctest -j16` (52 tests, ~9s).
- T1 tier-a: `python3 util/parity.py --old /tmp/opencode/om-ref/A0/openMalaria --new build-p5/openMalaria --recipe /tmp/opencode/om-ref/A0/recipe.json --tiers a -j16` (53). Tier-b: `--tiers b --seeds 1-10`.
- bench: `sh /tmp/opencode/om-bench/bench.sh /home/razlock/git/fork/build-p5/openMalaria TAG` (ABSOLUTE path!).
- perf: hardware events BLOCKED (perf_event_paranoid=2, no CAP_PERFMON) — use `perf record -e task-clock -F 3997` (software clock, no -g; dwarf callgraphs also blocked). Profile of p4b binary at /tmp/opencode/p5-perf/p5-sw.data: hot = dgemm 18% (CalcNvO GSL, phase-8 flip kills), getPEpisode 7.6%, probTransmissionToMosquito_falciparum 6.5%, update_descriptive+determineDensities ~9% (Infection loop).

## DONE
- **PERF FIX (step 1)**: Infection per-kind cold state moved to heap payloads (commit this): `EmpiricalData`/`MolineauxData`/`PennyData` behind `empirical`/`molineaux`/`penny` unique_ptrs; hot object ~500B→~100B for the benchmark kind (Descriptive, example@2000 uses default Descriptive model — profiled). Payloads allocated in ctorBody_X + read_X ONLY (no RNG draws there → m_rng stream untouched; wire byte-order UNCHANGED — checkpoint fns route through payloads, so 12/12a + all two-run box tests stayed green with NO golden touch). Dead `setPatentGrowthRateMultiplier` deleted (never called in old code either). `Kn` (PkPd IC50^slope cache, used by LSTMDrugType/Conversion via `inf->Kn`) → lazy `std::unique_ptr<map>` behind `Kn()` accessor (was 48B inline in every infection; NEVER serialised). unittest PennyInfectionSuite reads `infection->penny->threshold_*` (friend, same commit). Gates: ctest 52/52, unittest 91/91, tier-a 53/53, bench **12.48s mad .010** (−2.3% vs 12.77, floor 1.0% → KEEP), RSS 39.4→35.4MB. LOC model 18,599→18,614 (+15, payload decls — perf trade, honest ledger row). Residual vs A0: +1.2%.

## Queue
2. CYCLE BREAK #1: interventions bases out of model/Transmission/PerHost.h — census first (what forces the include), include-hygiene before any seam type.
3. CYCLE BREAK #2: dynamic_cast pokes in model/ → declared API (list sites in commit msg).
4. CLINICAL SEAM: formalize clinical→interventions (§5 one seam), document for phase 7.
5. FINAL GATE: T1 a+b 63/63, ctest, unittest, bench p5-final, cloc, census statics refresh if statics moved, PERF.md row 5.

## Learned / traps (inherited + new)
- m_rng(ctx.master_RNG) CONSTRUCTS-BY-DRAWING — only where old code did (4b h10 trap).
- ext lifetime only after ensure_extensions; no header-inline deref of ctx.ext; TwoRunsSuite = ext tripwire, rerun whenever RunExtensions changes.
- INHERITED-SCOPE lookup: qualify statics per body after merges; convert from `git show HEAD:file` anchors.
- Error strings byte-identical; `// parity:` comments travel with code.
- ESCaseManagementSuite.h + AnophelesModelSuite.h remain excluded+stale (pre-existing); deferred: setHS-overwritten clinical statics → ctx.ext.
- Untracked PROMPT.md/email = user files, leave.
