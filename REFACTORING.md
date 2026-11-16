# OpenMalaria Refactoring — qwen-refactor branch

Independent experiment, branched from `main` at `cfc3c8e0` (schema-50.0).
Goal: a simple, modern, KISS code base — plain data structures, free functions
that process them (`update(model.humans, model.transmission)`), no circular
dependencies, no global variables. C++20 is already enabled
(`CMAKE_CXX_STANDARD 20`).

This is NOT a port of any other refactor branch. Everything below was derived
from reading this tree.

---

## 1. State of the code today

~32,300 lines of C++ in `model/` (150 files), plus a build-time-generated XML
binding (`schema/` → `scnXml`, produced by `xsd` from `scenario_50.xsd`).

### The good — keep it

- **The test harness is excellent and is the contract.** `test/run.py` runs
  ~50 scenarios and compares survey output, continuous output AND resume-from-
  checkpoint files byte-for-byte against `test/expected/`. `unittest/` has
  Boost suites with golden `.csv`/`.yaml`/`.txt` files. Every phase below ends
  with these green, unchanged.
- **No classic singletons survive.** `Parameters` is already a value object
  built in `main.cpp` and passed by const-ref. Per-human RNG (`LocalRng`) is
  already properly derived from one master seed and threaded through calls.
- **The checkpoint mechanism** (`operator&` symmetric istream/ostream streams,
  `model/util/checkpoint.h`) is compact and format-stable. The byte layout, not
  the C++ plumbing, is the thing we must preserve.
- **Directory structure is honest**: Clinical, Host/(WithinHost), Transmission,
  PkPd, interventions, mon, util. Module boundaries exist; they are just
  breached everywhere.
- XML validation is generated from XSD. Never hand-roll XML parsing; keep
  `scnXml` as the single parsing layer and push it *behind* the model instead
  of removing it.

### The bad

1. **Time is a global.** `sim::now()/ts0()/ts1()/nowOrTs0()` are ambient
   reads: ~169 current-time reads in 69 of 150 files (`model/sim.h` statics).
   Functions that are conceptually pure take no time argument and read the
   clock from the ether. This is the deepest coupling in the code base.
2. **~25 "all-static class" pseudo-singletons.** A class of static functions
   with private static mutable data, `init()`-ed once from `main.cpp`, read
   from anywhere: `InterventionManager`, `Genotypes`, `PerHostAnophParams`
   (+`PerHost::relAvailAge`), `WithinHost::Treatments`, diagnostics,
   `AgeStructure`, `NeonatalMortality`, `InfectionIncidenceModel`, the
   Clinical/EventScheduler static clusters, `ITN/IRS/GVI/VaccineComponent`
   registries, PkPd drug registries, `mon::runtime`, `ModelOptions`,
   `CommandLine`, `master_RNG`, `StreamValidator`. ~95 mutable static members
   plus ~45 namespace-scope globals. Init order is an undocumented hand-tuned
   DAG of numbered comments in `main.cpp:152-194`.
3. **The Host ↔ Transmission ↔ Interventions cycle** (the reason this branch
   exists):
   - `Host::Human` embeds `Transmission::PerHost` (`Human.h:101`), so every
     human drags in the transmission headers.
   - `PerHost::update(Host::Human&)` calls back into the human (cohort
     removal) and is invoked from *inside* `VectorModel::calculateEIR`
     (`VectorModel.cpp:375`) — the net-hole tick is hidden in an EIR
     calculation.
   - `PerHost::hasActiveInterv` and every intervention effect getter route
     through the static `InterventionManager` + per-class `componentsByIndex`
     registries (`PerHost.cpp:137`, `ITN.cpp` 20+ sites).
   - `InterventionManager::init/deploy` takes `TransmissionModel&`;
     `TransmissionModel::updateKappa` iterates `Host::Human` and reaches into
     `withinHostModel`. `VectorModel` declares `friend class PerHost`.
   - Header cycle: `Human.h → PerHost.h → interventions/Interfaces.h → …` and
     20 headers `#include <schema/…>`, so the generated XML code leaks into
     almost every TU.
4. **The interventions subsystem is the worst area** (detailed in §5):
   duplicated age/entomo/coverage parsing (3-4 copies), 8 `TimedDeployment`
   subclasses that are one data struct + an executor, `Deployments.h` and
   `HumanInterventionComponents.h` are .cpp files in .h clothing (included by
   exactly one TU), deployment constructors mutate the vector model, four
   hand-rolled parallel component registries, ITN and IRS each re-implement
   the same survival/attractiveness math, `dynamic_cast` chains as poor-man's
   pattern matching, one **unchecked** dereference (`ITN.cpp:655`).
5. **Silent checkpoint gaps.** `AnophelesModel::operator&` does not serialize
   baited traps, non-human-host deployments, sugar-feeding or the
   string-keyed `reduce*` maps — resume-after-checkpoint silently loses
   vector interventions. `removeFromCohort` (`Human.cpp:131`) misses a cache
   invalidation that the ITN disposal path relies on.
6. **Dead weight**: `#if 0` blocks (R_0, life-cycle model), commented-out enum
   entries in `ModelOptions.h`, comments referencing files that no longer
   exist (`Interventions.cpp`), dead virtuals on `VectorModel` (no subclass),
   commented-out "cleanup" code in `transmission.h:359-365`.
7. **Old-C++ habits**: raw `new` in factory functions returning owning raw
   pointers, out-params where a return value fits, `map<string,…>` state,
   header monoliths, `cout` stringly-typed error assembly.

---

## 2. The hard constraints (what "byte-identical" really pins down)

1. **RNG draw order is part of the ABI.** Results are determined by the exact
   sequence of random draws: per human per time step, per species loop order,
   the inline copula/Bernoulli draws during deployment filtering, the master
   RNG draws in the mosquito model, even the "one extra draw for imported
   infection" quirk in `InfectionIncidenceModel`. Any restructure must keep
   the *interleaving* of draws identical — e.g. the net-hole Poisson draws
   happen per-human inside the human update loop (`calculateEIR`), not in a
   population pass. Moving a sampling call to a different loop level breaks
   the contract, even if the code looks nicer.
2. **Floating-point operation order is part of the ABI.** e.g. intervention
   effect factors are multiplied in `activeComponents` insertion order
   (`PerHost.cpp:104-131`), including not-yet-deployed components that return
   1.0. Keep fold order; do not "optimize" by skipping undeployed components
   or reordering factors.
3. **Checkpoint streams must round-trip.** `test/run.py` saves a checkpoint
   mid-run and resumes; survey + cts output after resume must match. No test
   loads a checkpoint from a *different* binary version (`test/old/` contains
   scenario XMLs only), so the stream layout may change — but only in commits
   where we deliberately change it, read/write stay perfectly symmetric, and
   the resumed outputs stay byte-identical.
4. **Output formatting** (survey `.sum`, cts files, scenario.sum) is compared
   verbatim. `mon/` formatting code is formatting, not model logic — touch
   with care, not at all where avoidable.
5. Schema version 50 XML in, identical numbers out — including all the
   "backwards-compatible bugfix" option paths in `ModelOptions`. Those
   branches are load-bearing for reproduction; simplify *where the option is
   read*, never collapse an option's behavior.

Safety net for every step: `ctest` in the build dir (51 tests: the curated
scenario suite with checkpoint save/resume + the Boost unit suites, all
byte-compared) — all green at baseline. NOTE: running `run.py` bare executes
53 scenarios, but `12a` and `MonitoringAll` are not in the suite lists in
`test/CMakeLists.txt` and have no expected files; use `ctest`, not bare
`run.py`, as the gate.

---

## 3. Target architecture

One owner of everything, plain data, free functions. `main.cpp` becomes:

```cpp
Model model = loadScenario(path);   // parse XML → Model (all reading here)
warmup(model);
while (model.clock.now() < model.end)
    step(model);                    // deploy → update → kappa → monitor
```

```cpp
struct Model {                    // model/Model.h — data + owned sub-models only
    SimClock      clock;          // replaces sim::now/ts0/ts1 globals (converters stay free funcs)
    Options       options;        // value struct; replaces ModelOptions/CommandLine statics
    MasterRng     master;
    Genotypes     genotypes;
    AgeStructure  demography;
    Population    humans;
    Transmission  transmission;   // one struct, mode enum + variant-like storage
    Interventions interventions;  // components vector + deployment specs (owned)
    Monitoring    mon;
};
```

Principles:

- **Free functions are the default; a virtual class must earn its keep.**
  Decided-at-startup → plain data struct + free function dispatch (a switch
  is honest where 8 TimedDeployment subclasses are not). Changes-at-runtime
  → virtual class (AnophelesModel vs SimpleMPDAnophelesModel, EmergenceRate
  estimators, WithinHost falciparum/vivax, Infection models, clinical
  decision trees). Small types with invariants (RNG, checkpoint streams)
  keep member functions: the invariant is the API.

- **Dependencies point one way**: `Transmission` never names `Host::Human`
  (it takes a `PerHostData&` + small `Cohorts&` view); `Host` never names
  `Transmission` types beyond a forward-declared interface passed to free
  functions; `interventions` read functions are the only code that sees
  `scnXml` (per-subsystem `readX(...)` free functions in `.cpp` files —
  no generated header in any `model/` header).
- **Polymorphism kept only where it earns its keep, and frozen at init**
  (user's call, decided once at scenario load):
  - `AnophelesModel` + `SimpleMPDAnophelesModel` → keep one virtual base,
    mark overrides `final`.
  - `EmergenceRateEstimator` (2 impls) → same.
  - `TransmissionModel` (NonVector | Vector): replace the base class with a
    mode enum + members; the two `calculateEIR` variants become two free
    functions chosen by the enum. Removes `friend class PerHost`, the dead
    virtuals, and 11 `dynamic_cast<VectorModel*>` sites (they become a
    `isVector(model.transmission)` check with the compiler enforcing it).
  - `PerHostInterventionData` (ITN/IRS/GVI): replace with one tagged struct
    per deployed component holding a `const ComponentParams&` (params live
    as long as the run — registries die) and a small state union/variant;
    effect factors become free functions. Serialization keeps writing the
    id first, so the per-type stream layout is untouched.
  - Clinical CM5Day/DecisionTree and WithinHost/Infection variant families:
    virtual interfaces are fine there (many implementations, genuine
    dispatch); only their *static parameter state* moves into param structs.
- **8 `TimedDeployment` subclasses collapse to data**: one
  `Deployment { Time date; AgeRange; EntomoRange; CohortFilter; Coverage(+copula);
  Action action; }` with `Action = std::variant<DeployHumanIntervention,
  ChangeHealthSystem, ChangeEIR, UninfectVectors, DeployVectorIntervention,
  DeployTrap, …>` and one `deploy(model, deployment)` free function. Repeat
  schedules stay an expanded multimap (it works and is checkpointed-by-replay).
- **Every init-time "static class" becomes a parameter struct** filled by a
  `readX(const scnXml::&, Options&)` function and stored in `Model` (or in the
  subsystem that owns it). Read paths take `const Params&`.
- **No `friend`, no `private`/`protected` data.** Every class is plain
  data with public members; access control is not doing useful work in a
  single-owner simulation model. RNG and checkpoint types keep member
  *functions* (the API), but their state is public. Enforced tree-wide
  (no `friend`/`private:` left in `model/**`).

- Checkpoint plumbing: keep `operator&` streams (the byte format is the
  contract), but each class gets ONE `template<class S> void ser(S&)` used in
  both directions — this deletes the hand-mirrored read/write field lists
  that are already drifting.

---

## 4. Phases (each one ends with byte-identical run.py + unittest)

**Phase 0 — baseline** (done): branch `qwen-refactor` from `main`; full build;
`test/run.py` + `unittest` green; baseline recorded here.

**Phase 1 — break the PerHost knot.** The user-visible coupling, first.
- `PerHost` loses its statics: `PerHostAnophParams::params` moves into the
  transmission model (species params it already owns); `relAvailAge` moves
  into a transmission-level `AvailabilityParams` value.
- `PerHost` becomes plain data (`struct PerHostData` — same fields, same
  checkpoint order) + free functions `Transmission::…`. It stops including
  `Host/Human.h`: the net tick takes a small `Cohorts&`-style view instead of
  `Host::Human&`; the tick call site moves from inside `VectorModel::
  calculateEIR` to an explicit hook in the per-human update — placed where
  RNG order is provably unchanged (first statement of `calculateEIR` ⇒ first
  thing before `getEIR` in `Host::update`).
- Intervention effects: `PerHostData` stores `const HumanVectorIntervention
  Component&` next to each deployed state; `hasActiveInterv` and effect
  getters never touch `InterventionManager`.
- End state: `Human.h` no longer includes PerHost internals; transmission
  compiles without `InterventionManager.cpp` statics.

**Phase 2 — `Model` struct + `step(model)`.** De-static
`InterventionManager`, `mon` runtime, `master_RNG`, `StreamValidator`,
`Genotypes`; they become `Model` members. `main.cpp`'s numbered init
comment-block becomes actual structure (the constructor order of `Model`).
`Host::update(human, transmission)` → `update(model)`.

**Phase 3 — kill the time global.** `sim::now/ts0/ts1/intervTime/intervDate`
reads become a `const SimClock&` (or `model.clock`) parameter. ~170 sites,
mechanical, per-directory commits. `sim::fromDays/inYears/oneTS`-style pure
converters stay exactly as they are (free, inline, const).

**Phase 4 — interventions rewrite** (the "untangle the mess", full rewrite):
- `Deployment` data model per §3; delete the 8-class hierarchy.
- One shared `readAgeRange/readEntomoFilter/readCoverage(+copula)` used by
  every deployment kind (currently 3-4 divergent copies — one of which sums
  the *raw* availability where the others sum the scaled one, and one has a
  mislabeled error string; normalize all to the behavior the scenarios pin).
- Shared `factors::` (survival/attractiveness blending) used by ITN and IRS
  instead of two copies; ITN net state machine stays semantically identical.
- Vector-side interventions: one checkpointed `VectorInterventionState`
  (replaces five string-keyed `reduce*` maps mutated from deployment
  constructors) — fixes the silent checkpoint gaps as part of the rewrite.
- Real `.cpp` files replace `Deployments.h` / `HumanInterventionComponents.h`.
- Delete the four `componentsByIndex` registries, the unchecked `dynamic_cast`
  at `ITN.cpp:655`, the stale-cohort-cache path.

**Phase 5 — XML interface separation.** All `scnXml` includes move into
`read*.cpp` translation units (one per subsystem: scenario root, entomology,
interventions, health system, monitoring, drugs, vivax). `ModelOptions` becomes
a value struct parsed once (bitset copy), with per-subsystem `Options` views;
the 59 scattered `ModelOptions::option()` reads become init-time caches (the
header already recommends this). 0 generated-header includes in `model/**.h`.

**Phase 6 — remaining static clusters.** Clinical (ClinicalModel,
EventScheduler ×14 statics, CM5DayCommon, ImmediateOutcomes,
`decision_library`), WithinHost (Treatments, diagnostics, Infection statics,
WHVivax), PkPd registries, InfectionIncidenceModel, Pathogenesis,
AgeStructure, NeonatalMortality → param structs owned by `Model`/subsystem,
filled by `read*` functions.

**Phase 7 — sweep.** Delete `#if 0` blocks, commented-out enum entries and
"cleanup" corpses, stale comments, dead virtuals; raw `new`/`unique_ptr`
factories → `make_unique`; out-params → returns where free; keep names,
keep includes minimal. Trivial C++20-isms only where they delete characters
(`span`, designated initializers); no style churn for its own sake.

---

## 5. Progress log

Phase 1 — break the PerHost knot — DONE (52/52 byte-identical throughout):

- 1.1 `b075639e` ITN/IRS/GVI components pass themselves (`*this`) into
  `makeHumanPart`; `componentsByIndex` static registries gone; unchecked
  `dynamic_cast` in `HumanITN::redeploy` gone (component reference stored).
- 1.2 `32ea03e6` `PerHostAnophParams::params` / `entoAvailabilityPercentiles`
  statics → `Transmission::VectorSpeciesParams` value member of
  `TransmissionModel`, threaded through `PerHost::initialise`, `Human`
  (via `Population`), `AnophelesModel::initAvailability/deployVectorTrap`
  and the entomo filters in `Deployments.h`. `main.cpp` builds the
  transmission model before the population (neither ctor draws RNG, so
  draw order is untouched).
- 1.3 `ecf10b64` `PerHost::relAvailAge` static curve → public
  `TransmissionModel::ageAvailability` (+ `setAgeAvailability` deriving
  `adultAge`); `PerHost` holds a non-owning pointer.
- 1.4 `0b3aa03d` net-hole tick: `VectorModel::calculateEIR` first line →
  `TransmissionModel::getEIR` first line via virtual `startHumanUpdate`
  hook (base no-op, vector overrides). RNG sequence position unchanged.
- 1.5 `6ee0a202` `getEIR`/`updateKappa` moved to new
  `TransmissionModel.cpp`; `TransmissionModel.h` no longer includes
  `Population.h` / WithinHost headers.

Phase 2 — de-static — in progress:

- 2.1 `fedd40fd` `util::master_RNG` → `MasterRng` local of `main()`,
  threaded by reference into transmission/population/checkpoint paths.
- 2.2a `5131e025` `PerHost::hasActiveInterv` answers from the deployed
  part itself (new virtual `componentType`), no registry lookup.
- 2.2b `b83f2d29` `InterventionManager` class-of-statics → plain
  `struct InterventionManager` (deployment state) owned by `main` +
  free `interventions::init/deploy/checkpoint/loadFromCheckpoint`.
  Component registry still file-static (Phase 4 replaces it).
- 2.3 `d4ef2ca5` `struct Model { master_rng, transmission, humans,
  interventions }` (model/Model.{h,cpp}) constructed in dependency order;
  `main()` owns one value, `run()`/checkpoint read/write take `Model&`.

Phase 4 (deployments) — first slice done:

- 4a `7278ce45` the `TimedDeployment` class hierarchy (Dummy/ChangeHS/
  ChangeEIR/UninfectVectors/Human/CumulativeHuman/VectorPop/Trap/ReduceNhh/
  AddNhh) and the `ContinuousHumanDeployment` class are decided at
  XML-read time → `struct TimedDeployment { date; variant<…actions> }` +
  free `apply/checkDeployDate/initReduceNhh/readAddNhh/filterAndDeploy/
  print_details`. `Deployments.h` 712 → 153 lines.
- 4b `9db62d02` one shared `inAvailabilityRange` for timed/cumulative/
  continuous (raw-vs-scaled is an explicit parameter, behaviour kept).

- 4c `56fd8cf1` decision-tree equivalence machinery deleted: `operator==`
  existed only for `save_decision`'s duplicate elision; both gone (15
  hand-written comparisons, ~165 lines). Pointer identity of trees was
  never observable (not checkpointed, exec is structural).

Phase 6 — remaining statics, one cluster at a time:

- 6a `af8abdb5` `NeonatalMortality` statics → `Model::neonatal` member;
  `ClinicalModel::update` / `Host::update` take it by reference.
- 6b `4ee4e68d` component registry de-static: `struct
  InterventionRegistry { identifiers; components }` is now a member of
  `InterventionManager` (built by `init()`, immutable afterwards). The
  free `getComponent`/`getComponentId` take it by const reference, and
  it is threaded explicitly to every construction site that resolves
  component ids: `HumanIntervention`/`TriggeredDeployments`/`Treatments`
  ctors, the whole `CMDecisionTree::create` chain (`CMDTDeploy`,
  `CMDTCohort`), `ScreenComponent`, `DecisionTreeComponent`,
  `ClinicalModel::setHS` → `DecisionTree5Day`/`ESCaseManagement`
  `setHealthSystem` → `WHInterface::addTreatment`, `mon::initCohorts`,
  and the checkpoint read path (`PerHost::load(istream, registry)`
  replaces the templated `operator&`, which stays ostream-only).
  File-statics `identifierMap`/`humanComponents` deleted.

Deliberately left file-static:

- `removeAtIds[SubPopRemove::NUM]` (HumanComponents.h) — written only
  during `init()`, read from deep inside `Human.cpp`; an init-only table,
  threading it through `Human::update` buys nothing yet.
- `decision_library` (CMDecisionTree.cpp) — arena owning decision-tree
  nodes; it *grows at runtime* (changeHS replays create new trees and
  old trees may still be referenced), so moving it into `Model` means
  threading an arena reference through every `create()` call on top of
  the registry. The tree interpreter keeps its virtuals by the §3 rule.
- `ModelOptions` — set once from scenario/CLI before the simulation and
  read (never written) at ~60 sites afterwards; a startup constant, so
  threading it as a value is churn without any state to protect.

Contract-blocked behaviour notes (fixing these changes outputs; needs a
decision to regenerate expected files, same class as the §6 checkpoint
gaps):

- `Human::removeFromCohort` (ITN fade-out) erases `subPopExp` but leaves
  the `cohortSet` monitoring bit set; episodes additionally report with
  the cohort snapshot taken at episode creation. Both are observable in
  survey output, so the asymmetry stays byte-identical.

Remaining Phase-6 clusters (not started, deliberately): the mon and
WithinHost static clusters. The clinical health-system statics are done —
see Phase 6c/6d below.

Phase 5 (XML interface at the boundary) — first slice:

- 5a `ec45f786` XML parsing moved out of widely-included headers:
  `Parameters.h`, `util/ModelNameProvider.h`, `util/ModelOptions.h`,
  `Transmission/TransmissionModel.h`, `Host/ImportedInfections.h` and
  `Clinical/ESCaseManagement.h` no longer include any generated schema
  header. New `Parameters.cpp`, `util/ModelNameProvider.cpp`,
  `Host/ImportedInfections.cpp` hold the moved XML-reading bodies;
  `ModelOptions::init`/`Parameters` ctor now take `const scnXml::Model&`
  (forward-declared) at the boundary. `averageEIR` and
  `setAgeAvailability` are out-of-line in `TransmissionModel.cpp` (both
  `averageEIR` overloads kept — `transmission.h` uses the DailyValues
  one); `schema/pharmacology.h` was an unused include in
  `ESCaseManagement.h`.

Phase 2 (end state) — finished:

- `61f46b86` `Host::update(human, transmission, neonatal)` →
  `Host::update(human, model)`.
- `6ec2f3bc` `interventions::deploy` / `loadFromCheckpoint` take
  `Model&`. The run loop and checkpoint resume now hand the whole model
  to the few functions that need several members.

Phase 7 (sweep) — first pass:

- `ac77e16e` all three `#if 0` R_0 corpses deleted (the live
  `insertR_0Case` scenario error stays); seven commented-out
  `ModelOptions` enum entries with their doc comments deleted; fifteen
  polymorphic `unique_ptr<Base>(new Derived(…))` factory sites converted
  to `make_unique<Derived>`.

Phase 7 (sweep) — dead-code hunt (−595 lines, 52/52 throughout):

- `cab6d3b3` deleted `DescriptiveIPTInfection` (never listed in
  CMakeLists, never referenced); free functions/macros/methods with zero
  call sites anywhere: `util::addTo`, `RNG::gen_seed`, `OM_ASSERT_EQ(_TOL)`
  + `do_assertion(_tol)`, `addProphylacticEffects` stub,
  `samplePatentValue`, `numDrugTypes`, `argmax`, `initIterateScale`,
  `IPTattenuateAsexualDensity` virtual, `VivaxBrood` `hasHad*` getters
  (the members stay — they are checkpointed).
- `80b993c6` commented-out code blocks deleted (old treatment-probability
  formula, hypnozoite partial-clearance alternative, GSL spline
  interpolator, disabled Molineaux tests, debug life-cycle loop) and the
  three never-registered, fully-commented `VectorModel::cts*` callbacks.
- `960e2073` `3fe8d99a` unused assert macros and the stale LCEmergence
  comment; the hardcoded-on `USE_DERIVATIVE` root-solver `#else` branch
  and the never-defined `WHVivaxSamples` debug blocks deleted (alternate
  GSL solver methods stay listed as comments — they are tuning knobs).
- `e4d8fc69` `-Wall -Wextra` fallout: counters set but never read
  (`nImported`, `sumLogDens`, AddNonHumanHosts `instance`), unused
  `deploy()` aliases, eight `return std::move(copy)` → `return copy`
  (NRVO) in `DecayFunction`.

Deliberately kept: `ModelOptions` flags whose only remaining reference is
the XML name table (old scenario files may still set them — deleting the
enum would turn a silent no-op into a scenario error), `StreamValidator`
(option-gated checkpoint tooling used by every `operator&`), the
net/IRS-insecticide TODO comments inside the empty-state callbacks they
refer to.

Phase 8 (access control) — friends and private data gone (52/52 throughout):

- Three commits: interventions headers first, `util/` headers next
  (RNG/Xoshiro state, samplers, decay wrappers, `traced_exception`), then
  every remaining header plus the local command/composite classes in
  `CMDecisionTree.cpp`, `DecayFunction.cpp`, `WHVivax.cpp`. All
  `private:`/`protected:` labels became `public:`, all 34 `friend`
  declarations deleted (test-suite friends were pure access grants),
  `SimData` became a plain struct. Pure access change — layout, names,
  draw order untouched.

Phase 5b (XML boundary) — no model header outside the boundary includes
the schema (52/52 throughout):

- `04862052` `sampler.h` (the `createSampler<T>` template only served two
  samplers → two file-local builders in `sampler.cpp`, every error message
  preserved); `d645d342` ITN/IRS/GVI/Vaccine — the nested XSD spellings
  (`ITNDescription::AnophelesParamsType`) turned out to be typedefs of
  top-level classes (`scnXml::AnophelesParams1` etc.), so forward
  declarations suffice; `b5ade779` `InterventionManager`/`Deployments`/
  `Interfaces`/`LSTMMedicate` + CMake `-DXSD_CXX11` (generated headers
  defined it themselves before the XSD prelude; making it global lets
  `containers.hxx` be included in any order); `b22354a7` deleted
  `HumanInterventionComponents.h` (definitions-only header included by
  exactly one source — inlined into `InterventionManager.cpp`) and cleaned
  the clinical headers. Removing the old includes also removed crutches:
  several transmission/clinical headers had relied on transitive schema
  includes and now forward-declare or include what they parse themselves.

Status of the remaining plan (what is *not* done yet):

- **XML boundary (Phase 5)**: DONE. 5a moved core-header parsing to
  `.cpp` files; 5b removed the schema include from every remaining target
  header: `sampler.h` (XML construction → two builders in `sampler.cpp`),
  `ITN`/`IRS`/`GVI`/`Vaccine` (nested XSD spellings like
  `ITNDescription::AnophelesParamsType` are just typedefs of top-level
  classes — forward-declare the real name instead), `InterventionManager`,
  `Deployments` (plus `initReduceNhh`/`readAddNhh` became file-local in
  `InterventionManager.cpp`), `Interfaces`, `LSTMMedicate`, and
  `CMDecisionTree`/clinical headers. `HumanInterventionComponents.h` —
  definitions only, included by exactly one source — was inlined into
  `InterventionManager.cpp` and deleted. The three entomology-factory
  headers were the last tail: `transmission.cpp` (`20a63ea5`) took the
  inline factories, `NonVectorModel.cpp` (`6d2fc059`) its ctor/`changeEIR`
  bodies, and `3a4130e3` replaced `VectorModel`/`AnophelesModel`'s
  nested-sequence signatures with pre-parsed `VectorPopIntervDesc`/
  `VectorTrapDesc` structs built at the `InterventionManager` boundary.
  Schema headers are now included only by `DocumentLoader.h`, `XMLChecker.h`
  and the `.cpp` files that actually parse XML — **no model header includes
  the schema anymore.** Global `-DXSD_CXX11` in CMake lets
  `xsd/cxx/tree/containers.hxx` be included in any order (the generated
  headers always defined it themselves before their prelude).
- **Interventions rewrite (Phase 4)**: 4a (deployment class hierarchy →
  variant data + free functions), 4b (one shared entomo filter) and 4c
  (one shared `factors::` survival/attractiveness model for ITN and IRS —
  IRS sets the hole factors to 0 and evaluates at holeIndex=0, which is
  bit-exact, plus a `clamp` flag for ITN's out-of-bounds handling) done.
   4d collapsed the vector-intervention components: `PerHostInterventionData`
   is one Kind-tagged struct with per-kind free functions, and ITN/IRS/GVI
   components are one concrete `HumanVectorInterventionComponent` (Kind +
   public data, parse functions at the boundary). Left in this family: the
   six trivial IM-local component classes and `VaccineComponent` — cold-path
   `deploy()` only, low value; not started.
   The silent checkpoint gap is CLOSED (see below): the five string-keyed
   `reduce*` maps, `nhhInstances`, `baitedTraps` and
   `probAdditionalDeathSugarFeedingIntervs` are now checkpointed.
- **Phase 3 (clock global)**: in progress — `SimClock`/`sim::now` etc.
  conversions landed for `mon/Continuous`, the whole human-intervention
  deploy chain, vaccine decay, the per-host intervention/effect chain
  (`model/interventions/` has zero clock-global reads) and the full
  transmission update chain incl. `AnophelesModel` and the `Human`/
  `Population` step (`b96bb10a`; `TransmissionModel.h`/`NonVectorModel.h`
  dynamic reads and `AnophelesModel.cpp` are clock-only). Remaining reads
  (~105, incl. `main.cpp`/`checkpoint.cpp`/`mon` themselves): the
  within-host/clinical models (mon/WithinHost cluster), the
  deferred cts-callback/init-time reads listed under `b96bb10a` below, and
  `main.cpp`/`sim` themselves (the owner — dies with the last slice).
- **Phase 6**: `NeonatalMortality` (6a), component registry (6b), clinical/HS
  statics (6c, `Clinical::ClinicalParams`), `InfantMortality` (6d),
  `InfectionIncidenceModel` (6e, also de-virtualised) and the PkPd
  registries (6f, `PkPd::PkPdState`) done; the mon/WithinHost and
  Pathogenesis static clusters remain.
- 59 `ModelOptions::option()` reads remain (Phase 5 wants init-time caches).

Net effect so far: no mutable process-wide state left in the
transmission/PerHost path; intervention deployment state, the component
registry, the master RNG, the neonatal and infant-mortality models and all
clinical/health-system parameters are caller-owned (`Model` members). The XML
boundary is complete — schema headers reach only the loader and the `.cpp`
files that parse. Cumulative branch diff vs `main`: +6577/−6876 across 144
files — the model tree is 31,226 lines.

---

Phase 4c (shared entomo factors) — one survival/attractiveness model for
ITN and IRS (−209 lines, 52/52):

- `factors::SurvivalFactor`/`RelativeAttractiveness` moved verbatim from
  `ITN.h`/`ITN.cpp` into `interventions/InterventionFactors.h`; the IRS
  nested copies were the same equations with hole terms absent, so IRS
  initialises the shared structs with `HF=IF=lHF=lIF=holeScaling=0` and
  evaluates at `holeIndex=0` — the extra terms are exact IEEE no-ops, so
  results are bit-identical — plus one `clamp` flag distinguishing ITN's
  out-of-bounds clamping from IRS's assert-only path.

---

Phase 3 (clock global) — slices (52/52, bit-identical):

- `f69d61c4` added `struct SimClock{now, intervTime, intervDate}` and made
  `mon/Continuous::update` take it (one caller). `6d0489f8` threaded a
  `SimTime now` (== `sim::nowOrTs1()`) through `HumanInterventionComponent::
  deploy` and the `HumanIntervention`/`TriggeredDeployments` wrappers, and a
  `SimClock` through the deploy loop / timed-visitor / `loadFromCheckpoint`.
  `6c3f6c69` gave `PerHumanVaccine` the clock; `af125d80` threaded the clock
  through `PerHostInterventionData` redeploy/update/effect virtuals, the
  `PerHost` aggregators and `makeHumanPart`. Every caller passes the same
  value the code used to read internally, so all outputs are unchanged; the
  interventions directory no longer reads a time global.
- `b96bb10a` transmission update chain: `SimClock` gained `ts1` (==
  `sim::nowOrTs1()`), `main::run` builds a loop-top clock (continuous
  reporting, deploys, survey `summarize`) and a post-`start_update`
  `update_clock` (`now == ts0`, `ts1 == end of step`), threaded through
  `getEIR`/`calculateEIR`/`startHumanUpdate`/`updateKappa`/`surveyEIR`/
  `summarize`/`vectorUpdate`, `AnophelesModel::advancePeriod`/`update`/
  `summarize(species)`/`getInitPartialEIR`/`getInterventionEIR`/
  `getEmergenceRate` (now with an `intervTime` param), the vector
  intervention deploys (`deployVectorPopInterv`/`deployVectorTrap` from the
  timed-action `SimClock`), `Host::update`, `Population::update(ts1)` and
  `PerHostInterventionData`'s ctor (`deployTime` param; the checkpoint-load
  ctor value is overwritten from the stream). Deferred, documented:
  `VectorModel::init2` (init-time) and the continuous-report callbacks
  (`ctsCb*`, `getLastVecStat`, `getResAvailability`) — their
  `mon::Continuous` callback signatures carry no clock.

---

Phase 4d — vector-intervention component hierarchy → data (52/52 each):

- `e25780b6` — `PerHostInterventionData`'s three-class hierarchy (`HumanITN`/
  `HumanIRS`/`HumanGVI`, 16 virtuals) collapsed into one Kind-tagged
  public-data struct (per-host deployment state: decay heterogeneity sample +
  ITN hole/net fields + shared initial insecticide) plus 27 per-kind free
  functions (`itn*`/`irs*`/`gvi*` InitNew/Redeploy/Update/4 effects/
  Checkpoint/Read, declared in PerHost.h, defined in the kind's .cpp,
  dispatched by switch in PerHost.cpp's anonymous namespace). Per-kind
  checkpoint field order untouched.
- `6701e72c` — `ITNComponent`/`IRSComponent`/`GVIComponent` merged into one
  concrete `Transmission::HumanVectorInterventionComponent` (kind tag + flat
  fields; per-species `SpeciesParams` holds the ITN/IRS response curves and
  the GVI constant effects; GVI's `decay` folded into `insecticideDecay`).
  `deploy`/`print_details` implemented once (kind switch for the monitoring
  measures); `sampleInitialInsecticide` shared (identical clamp order →
  bit-exact). Each header shrank to its `parseXComponent(id, xml, species)`
  factory declaration — the XML boundary builds the struct, validation and
  error strings preserved. −122 lines net.

---

Phase-5 tail (XML out of the last three headers, 52/52 each):

- `20a63ea5` — `transmission.h`'s ~370 lines of inline entomology factories
  moved to a new `transmission.cpp` (anonymous-namespace helpers; the header
  now only declares `createTransmissionModel`); `checkpoint.h` dropped its
  now-unneeded schema-reachable include.
- `6d2fc059` — `NonVectorModel`'s ctor + `changeEIRIntervention` bodies moved
  to a new `NonVectorModel.cpp`; the ctor's unused `Entomology` parameter
  dropped; header keeps a `scnXml::NonVector` forward declaration.
- `3a4130e3` — the last nested-type signatures: `VectorModel::initVectorInterv`/
  `initVectorTrap` (and `AnophelesModel`'s consumers) now take pre-parsed
  `Anopheles::VectorPopIntervDesc`/`VectorTrapDesc` (plain data +
  `unique_ptr<DecayFunction>`, built by parse helpers in
  `InterventionManager.cpp`, validation and error strings preserved; trap
  name becomes `std::optional<std::string>`); `SimpleDecayingValue::set`
  takes the built decay function instead of XML. The two `init*` methods lost
  their never-overridden `virtual`. **No model header includes the schema now**
  (only `DocumentLoader.h`/`XMLChecker.h`, the boundary itself).

---

Checkpoint-gap fix (`15b30e78`, 52/52): the byte-identical checkpoint FORMAT
restriction was lifted (the output goldens remain the contract).
`AnophelesModel`'s stream now carries what deployments create but never
survived a resume: the five string-keyed `reduce*` NHH maps, `nhhInstances`,
`baitedTraps`, `probAdditionalDeathSugarFeedingIntervs`. Pattern: decay
parameters are rebuilt from the scenario XML before the checkpoint is read,
so only deployment state streams — reduce* maps overlay their existing
entries (sizes verified), traps/NHH rebuild in full with `validateListSize`
guards; `SimpleDecayingValue::operator&` and
`DecayFunction::checkpointSample` do the rest. Boxtests exercise the pair on
populated maps (NHH/reduce state exists at init; a misaligned stream
corrupts every later field and fails the golden comparison). Mid-run
deployment restore is symmetric by construction but exercised only by
inspection — the harness checkpoints at the init boundary only.

---

De-polymorphisation sweep (virtuals → Kind tags + switch; 52/52 each,
bit-identical by construction — same arithmetic, same RNG order):

- `380dfbc6` `Sampler`: the abstract base had exactly two consumer slots
  (`DecayFunction::het`, `PerHost::entoAvailability`). `LognormalSampler` and
  `GammaSampler` became plain final structs (kept their `unique_ptr`
  `fromMeanCV`/`fromMeanVariance` builders — every concrete caller stores
  `unique_ptr<LognormalSampler>` anyway); the base was replaced by a final
  `Sampler` wrapping `std::variant<Lognormal, Gamma>` with four
  `std::visit` members; gamma keeps a `sample(NormalSample)` that throws the
  same message the old base did. `DecayFunction::het` — always lognormal —
  became a value member (clone() gone).
- `66439abe` `DecayFunction`: 12 file-private subclasses (constant…emax,
  incl. the operator/Emax recursion that made `clone` virtual) collapsed into
  the one class: `Kind` enum, the union of plain parameters
  (`invL,k,constOverLambda,Emax,IC50,slope,initialConcentration`), child
  pointers `f1,f2`, deep-copy ctor, `compute`/`sampleAgeOfDecay`/
  `makeHetSample` switches (−102 lines). Checkpoint format untouched: only
  `hetFactor` was ever streamed; children re-cloned via `hetSample`.
- `86c19137` decision trees: 15 node classes with `virtual exec` became one
  public-data `CMDecisionTree` (`Kind` + all fields public; the Deploy node
  holds `optional<HumanIntervention>` instead of privately inheriting it)
  with `exec` a switch and the per-kind XML builders as file-local free
  functions. `save_decision`/`decision_library` unchanged.
- `a97f2566` `AgeGroupInterpolation`: only two implementations were live
  (linear, none) plus a throwing dummy; the whole hierarchy collapsed into
  `AgeGroupInterpolator` itself — `Kind{None,PiecewiseConstant,PiecewiseLinear}`
  + one `map<double,double>`, no heap, no dummy singleton.
- `26a473a4` `LSTMDrug`: `LSTMDrugOneComp`/`ThreeComp`/`Conversion` (18
  virtuals) → one Kind-tagged struct — flat fields for all three PK models,
  the virtual-method bodies as `lstmdrug::oneComp/threeComp/conversion` free
  functions in the same .cpp files (GSL workspaces and `Params_*` structs now
  file-local), dispatched by switch in `LSTMDrug.cpp`; per-model construction
  (sample order preserved) moved into `LSTMDrugType::createInstance`; the
  three-constructor ThreeComp cache went `mutable` so `drugFactor` stays
  const. Per-kind checkpoint field order byte-identical (−209 lines).

Remaining virtual clusters are the by-design strategy families (WH models,
vector/transmission/pathogenesis models) plus a handful of small component
classes (`InterventionManager.cpp`'s six trivial `HumanInterventionComponent`
kinds and `VaccineComponent` — cold-path deploy only, not worth a Kind tag).

---

Phase 6c — clinical/HS statics → one plain `Clinical::ClinicalParams`
struct (52/52 byte-identical, net −2 lines):

- `46013e1f` — new `Clinical/ClinicalParams.h`: one public-data struct with
  all former clinical statics — kind/option flags, shared HS data
  (`healthSystemMemory`, CFR + sequelae interpolators, `oddsRatioThreshold`),
  the `CM5DayCommon` case-management arrays and severe treatment,
  `ImmediateOutcomes` cure rates + diagnostic flag, the two
  `DecisionTree5Day` tree roots, and the whole `ClinicalEventScheduler`
  parameter set (durations, immunity curve, NMF log-odds/interpolators, its
  two tree roots). The `CaseType` enum moved to namespace scope. Ctor +
  `loadHealthSystem` + `getCommunityCFR` implemented in `ClinicalModel.cpp`.
- Static init/load became free functions: `ClinicalModel::init`/`setHS` →
  the `ClinicalParams` ctor + `loadHealthSystem` (dispatches to free
  `esInit`/`loadHSEventScheduler`/`loadHSImmediateOutcomes`/`loadHSDT5Day`
  in the kind files); the `ESCaseManagement` class → free
  `loadHSEventScheduler`/`esExecute`; `createClinicalModel(tSF, params)`
  reads kind flags from params.
- Ownership/threading: a value member of `Model` (declared between
  `transmission` and `humans`), constructed in `main` before the Model (the
  Human ctors call `createClinicalModel`); threaded `const ClinicalParams&`
  through `Population`/`Human` ctors, `ClinicalModel::update`/
  `doClinicalUpdate`/`isExistingCase`/`Host::summarize` and
  `Episode::update` (which now takes `hsMemory`); mutated only via
  `model.clinical_params.loadHealthSystem` (`main` after
  `interventions::init`, and the `ChangeHS` timed action). Not
  checkpointed — rebuilt from XML, exactly like the old statics.
- `CMDecisionTree::create` gained a trailing `SimTime hsMemory` (validates
  `<uncomplicated>` nodes; replaces the `ClinicalModel::hsMemory()` static
  read) threaded through the branching creator chain; `interventions::init`
  takes it for `DecisionTreeComponent` (the old static happened to be set
  by then); `CM5DayCommon`/ES keep their per-model virtuals — that
  polymorphism earns its keep (§3).

---

Phase 6d — InfantMortality statics → struct (52/52 byte-identical):

- `74f2fb74` — the file-static counters (`infantDeaths`,
  `infantIntervalsAtRisk`, `nonMalariaMortality`) and the `InfantMortality`
  static-method class became one plain struct in a new
  `Clinical/InfantMortality.h` (inline members, mirroring
  `Host/NeonatalMortality.h`), a value member of `Model`. `init` moved from
  the `ClinicalParams` ctor to `main` (next to `neonatal.init`);
  `ClinicalModel::updateInfantDeaths` takes the struct; checkpointing is its
  `operator&` at the identical stream position; `mon::writeSurveyData` now
  takes the `allCause()` value (computed in `main`, threaded through the
  file-local writer chain) so `Monitoring.cpp` no longer includes any
  clinical header.

---

Phase 6e — InfectionIncidenceModel de-static + de-virtual (52/52):

- `c01305e2` — the ~14 file-static parameters/option flags and the
  four-class `virtual` hierarchy (base + HeterogeneityWorkaroundII/
  NegBinomMAII/LogNormalMAII, no subclass data members — exactly the
  header's `TODO(optimisation)`) became: a plain
  `Host::InfectionIncidenceParams` (kind tag + doubles + two option flags,
  `init()` in main, value member of `Model`, threaded to the `Human` ctor
  and update calls) and one per-host state struct (`m_pInfected`,
  `m_cumulativeEIRa` — now a by-value `Human` member, checkpoint stream
  unchanged); the per-kind bodies are a switch in the free
  `getAvailabilityFactor`/`getModelExpectedInfections` (draw order kept).
  `ctsNewInfections` (checkpointed at its old stream position, in `Model`)
  is written by `reportNumNewInfections` and read/reset by a lambda
  registered with `mon::Continuous` from `main` (callback map is
  name-keyed, so registration order is not observable).

Phase 6f — PkPd registries de-static (52/52):

- `ee0d6e4e` — `PkPd::LSTMDrugType`'s and `PkPd::LSTMTreatments`' file-static
  lists/names became plain `DrugLibrary{types,names,inUse}` and
  `TreatmentLibrary{schedules,scheduleNames,dosages,dosagesNames}` with free
  `init`/`findDrug`/`findSchedule`/`findDosages`/`createInstance`; bundled as
  move-only `PkPdState{drugs,treatments}` (`LSTMDrugType` is move-only; the
  vectors heap-allocate, so type-object addresses stay stable for
  `LSTMDrug`'s pointers across the move into `Model`). `PkPdState::init`
  keeps the `pharmacology present()` check and is called from
  `WHInterface::init` under the same `opt_common_whm` condition that used to
  gate `LSTMModel::init`, preserving conditional init; `LSTMModel` now takes
  a `PkPdState&` (pointer member, stream layout unchanged). The library ref
  threads through `createWithinHostModel`→`CommonWithinHost`, and through
  `ClinicalParams::loadHealthSystem` → `loadHSEventScheduler`/`loadHSDT5Day`
  → `CMDecisionTree::create` (builders + `TreatInfo` take a trailing
  `TreatmentLibrary&`), and `interventions::init` → the
  `TreatPKPDComponent`/`DecisionTreeComponent` ctors plus `apply()` for the
  ChangeHS reload. `clear()` statics dropped; the unit-test suites hold
  their own `PkPdState` (`PkPdSuiteSetup` returns one, `WHMock` takes one,
  `createHuman` passes a null-state sentinel since it builds descriptive
  within-host models).

---

## 6. Intervention subsystem map (for the Phase-4 rewrite)

```
reading:      InterventionManager::init  (if/else chain over 13 component kinds)
storage:      identifierMap text→ComponentId; humanComponents vector; 4 registries
scheduling:   ContinuousHumanDeployment (per-human scan) + TimedDeployment
              multimap<SimTime, unique_ptr<…>>; replayed after checkpoint load
target:       HumanIntervention = sorted list of ComponentIds (+3× copied ctor)
per-human:    cohorts: Human::subPopExp map + cohortSet bit-cache (+1 stale path)
              PerHostData::activeComponents (ITN/IRS/GVI per-species effects)
              PerHumanVaccine (bypasses PerHost entirely — fine)
vector-side:  AnophelesModel string-keyed decay maps + traps + NHH instances
checkpoint:   ids written per object; dynamic_cast factory on read; timed lists
              rebuilt from XML, changeHS/changeEIR replayed; vector-side NHH/
              trap/reduce*/sugar-feeding state now checkpointed in
              AnophelesModel
```

What ITN actually computes (nothing more): per (human, ITN component):
insecticide amount (lognormal, clamped), disposal time (Weibull), hole count
(Poisson per step) and hole index (holes + ripFactor·rips, Poisson), then
at effect-query time evaluates a survival/attractiveness curve at
(holeIndex, insecticide at now) and blends protected/unprotected by
`propUse·propActive`. IRS is the same without holes. GVI is one decay value.
All of this is ~50 lines of real math wearing 741 lines of scaffolding.

---

## 7. Ground rules for every commit

1. Build clean; `unittest` green; `test/run.py` byte-identical. No exceptions.
2. One mechanical change per commit; behavior-preserving by construction
   (if a sampling call moves, prove the interleaving is identical in the
   commit message).
3. When a checkpoint stream must change, say so in the commit message, keep
   read/write symmetric, and rely on run.py's resume test.
4. Comment the *why* at every deliberate simplification, mark known ceilings
   with `// ponytail:` comments.
5. Deletions > additions. Every phase should shrink the tree.
