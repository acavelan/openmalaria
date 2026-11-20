# Phase 7f HANDOFF — Infection slim-struct + per-kind free functions (branch `p7f-infection` off tag `p7e-green` = b9f7c22a)

User-approved readability refactor (user vetoed re-litigating the 4b merged
god-class shape): `Infection` becomes a SLIM SHARED STRUCT + per-kind state
payloads + FREE FUNCTIONS (iv/vac/ComponentConfig pattern). NAMESPACING ONLY:
formula bodies move verbatim, RNG draw points/order frozen, checkpoint wire
order per kind identical (bytes kept), error strings verbatim, `// parity:`
comments travel.

## Environment
- Gate build: `build-p7f/` = Release `-O3 -DNDEBUG -ffp-contract=off`
  `-DOM_BOXTEST_BINARY=ON`. `cmake --build build-p7f -j16`.
- unittest: run from `build-p7f/unittest/` (`./unittest` → ".OK!", 94 cases;
  fixture files live next to the binary).
- ctest: `cd build-p7f && ctest -j16` (52).
- T1: `python3 util/parity.py --old /tmp/opencode/om-ref/A0/openMalaria --new
  build-p7f/openMalaria --recipe /tmp/opencode/om-ref/A0/recipe.json --tiers a
  -j16` (53); `--tiers b --seeds 1-4` (4). parity.py resolves
  `scenario_current.xsd` via `<repo>/schema`, `<repo>/test`, then
  `<repo>/build/schema` — `build/` here is the OLD real build dir left in place
  on purpose (p7e gates ran against its schema copy; do not delete).
- Bench (final only): `cd /tmp/opencode/om-bench && ./bench.sh <abs-bin>
  /home/razlock/git/fork/util/example/example_scenario.xml` — taskset -c 5,
  1 discard + 5× %U CPU-s; never concurrent with other load.
- Baseline stamp (p7e-green, build-p7f): ctest 52/52, unittest .OK!, tier-a
  53/53. cloc code LOC: model 18,417; unittest 3,794; Infection/ dir 1,091
  (1,891 total lines).

## Design (locked)
- `struct Infection` — public data, no friends (7.2b precedent "payloads
  public, friends dead"). Shared fields: kind, ctx, m_startDate, m_density,
  m_cumulativeExposureJ, m_genotype, m_origin, s_latentP (static, shared),
  payloads descriptive/empirical/molineaux/penny/kn (public unique_ptrs).
- DescriptiveData{duration,notPrintedMDWarning} payload joins the phase-5 trio;
  Descriptive statics (meanLogParasiteCount, sigma0sq, xNuStar), Empirical
  _* statics and Penny immune_threshold_gamma/update_density_gamma move to
  file scope of their owning TUs (grep-verified: no external readers; all
  init_* re-run per scenario load → per-run re-init preserved; TwoRunsSuite is
  the tripwire).
- Kn() stays SHARED on the struct: callers are LSTMDrugType.cpp:59-63 and
  LSTMDrugConversion.cpp:174-186 (PkPd LSTMDrug path); PkPd is unavailable for
  Descriptive (treatPkPd_descriptive throws) but live for ALL FOUR Common
  kinds → multi-kind → shared. Lazy unique_ptr kept (phase-5 allocation fix).
- Naming: shared dispatch keeps `infection_` prefix (infection_new,
  infection_new_checkpointed, infection_update, infection_checkpoint,
  infection_init_<kind>, infection_ctor_<kind>, infection_duration_descriptive);
  per-kind hooks are verb-first free functions (update_density_<kind>,
  read_<kind>, checkpoint_<kind>, determine_densities_descriptive). Declarations
  in Infection.h. TU-only helpers (getInflatedDensity, sigma_noise,
  sampleSubPatentValue, getClonalSummation, getVariantSpecificSummation) become
  file-local statics in their kind TU.
- Trivial accessors (getDensity, cumulativeExposureJ, clearImmunity,
  bloodStage, isHrp2Deficient, genotype(), origin(), expired(), seqDensity(),
  Kn(), operator&) stay inline members of the struct — zero WH churn.
- Ctor draw point: shared ctor keeps NO RNG draw (nowOrTs0() is a pure getter);
  Descriptive duration draw moves into `case Kind::Descriptive:` of the
  infection_new switch (infection_ctor_descriptive →
  infection_duration_descriptive draw is that case's FIRST and only RNG action;
  the assert before it never draws). Checkpoint reader for Descriptive
  allocates its payload before reading duration+notPrintedMDWarning (wire
  order unchanged: those two last, after the 5 shared fields).
- Commit ladder: 1) class→struct public + friends deleted; 2) Descriptive;
  3) Dummy+Empirical; 4) Molineaux+Penny; 5) shared ctor/update/checkpoint
  dispatch flip + m_kind→kind + header shrink + suites; 6) final gate + PERF
  row + close. Every commit: build + ctest 52 + unittest 94 + tier-a 53;
  tier-b seeds 1-4 after Descriptive and Molineaux commits; full a+b at final.

## Status
- [x] setup: branch, build-p7f green (52/52, .OK!, 53/53), this HANDOFF
- [x] C1 struct flip 6d1f0189 (public struct, friend blocks + MolineauxInfectionSuite/
      PennyInfectionSuite/UnittestUtil friend decls dead)
- [x] C2 Descriptive 64149ea0 (payload + statics→file scope + free fns + call sites)
- [x] C3 Dummy+Empirical 48b4e469
- [x] C4 Molineaux+Penny b0f16140
- [x] C5 shared dispatch flip ddce7ebc (infection_new/update/checkpoint), kind rename
- [x] C6 final: T1 a+b 63/63, cloc table, bench, PERF row 7f, report
- [x] perf-fix 6f025ebb: bench found +2.3 % regression (heap DescriptiveData
      payload: 2nd malloc/infection + dependent load per tick in expired());
      fixed by embedding the 12-byte payload by value. Re-bench ≈12.17 vs
      ≈12.24 (≤1 % band); T1 re-run green after fix.

### Close-out (7f final)
- Final gate on 6f025ebb: ctest 52/52, unittest .OK! (94), T1 tiers a+b
  seeds 1-10 = 63/63 bit-identical.
- LOC (cloc code): model 18,417 → 18,421 (+4); Infection dir 1,091 → 1,095
  (+4); header 417 → 402 lines. Readability refactor — LOC-neutral by design.
- Bench (taskset c5, median user-s, interleaved): p7e ≈12.15-12.16,
  p7f ≈12.23-12.25 → +0.5 % residual, perf-attributable to inline/layout
  shift only (determine_densities_descriptive 4.33→4.65 %, update_descriptive
  4.00→4.32 %); accepted (<1 % band).

## Traps (inherited + new)
- -ffp-contract=off always in gated builds; NO reassociation (phase 8 owns it).
- m_rng draw-order rule: the Descriptive duration draw, SlotAllocatorSuite
  state-word goldens ("genotype draw + duration draw") and RngGoldens are the
  tripwires; TwoRunsSuite catches statics-reset mistakes.
- unittest must run from build-p7f/unittest/.
- Scenario load calls Infection::init(latentP) + exactly one init_<kind>
  (model/xml/scenario.cpp:97-141) — keep call order on the free-fn retype.
- Checkpoints: bytes kept this phase (12/12a goldens stay green) — format may
  break per user ruling but this phase buys nothing from it.
- ≤5 fix iterations per divergence, then revert+log.
- Untracked PROMPT.md/email/pgodata-p7d/ = user files, leave alone.
