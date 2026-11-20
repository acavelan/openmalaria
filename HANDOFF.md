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
- [ ] C1 struct flip (public struct, friend blocks + MolineauxInfectionSuite/
      PennyInfectionSuite/UnittestUtil friend decls dead)
- [ ] C2 Descriptive (payload + statics→file scope + free fns + call sites)
- [ ] C3 Dummy+Empirical
- [ ] C4 Molineaux+Penny
- [ ] C5 shared dispatch flip (infection_new/update/checkpoint), kind rename
- [ ] C6 final: T1 a+b 63/63, cloc table, bench median-of-3, PERF row 7f,
      report (field layout, statics verdict, Kn placement, ctor-draw proof,
      phase-8 handoff note)

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
