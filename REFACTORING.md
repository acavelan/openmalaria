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

Remaining big Phase-6 cluster (not started, deliberately — no
half-threaded state exists): the clinical health-system parameters
(`CM5DayCommon`, `ImmediateOutcomes`, `ESCaseManagement`,
`EventScheduler`, `ClinicalModel` statics). They are rewritten by
changeHS at runtime and read from per-human clinical code and from
static decision-tree paths; threading ≈30 statics through ≈40
signatures is a multi-day pass of its own.

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

- `475e13ce` deleted `DescriptiveIPTInfection` (never listed in
  CMakeLists, never referenced); free functions/macros/methods with zero
  call sites anywhere: `util::addTo`, `RNG::gen_seed`, `OM_ASSERT_EQ(_TOL)`
  + `do_assertion(_tol)`, `addProphylacticEffects` stub,
  `samplePatentValue`, `numDrugTypes`, `argmax`, `initIterateScale`,
  `IPTattenuateAsexualDensity` virtual, `VivaxBrood` `hasHad*` getters
  (the members stay — they are checkpointed).
- `7b31ce21` commented-out code blocks deleted (old treatment-probability
  formula, hypnozoite partial-clearance alternative, GSL spline
  interpolator, disabled Molineaux tests, debug life-cycle loop) and the
  three never-registered, fully-commented `VectorModel::cts*` callbacks.
- `92d0db3d` `62eba769` unused assert macros and the stale LCEmergence
  comment; the hardcoded-on `USE_DERIVATIVE` root-solver `#else` branch
  and the never-defined `WHVivaxSamples` debug blocks deleted (alternate
  GSL solver methods stay listed as comments — they are tuning knobs).
- `c5e606b2` `-Wall -Wextra` fallout: counters set but never read
  (`nImported`, `sumLogDens`, AddNonHumanHosts `instance`), unused
  `deploy()` aliases, eight `return std::move(copy)` → `return copy`
  (NRVO) in `DecayFunction`.

Deliberately kept: `ModelOptions` flags whose only remaining reference is
the XML name table (old scenario files may still set them — deleting the
enum would turn a silent no-op into a scenario error), `StreamValidator`
(option-gated checkpoint tooling used by every `operator&`), the
net/IRS-insecticide TODO comments inside the empty-state callbacks they
refer to.

Net effect so far: no mutable process-wide state left in the
transmission/PerHost path; intervention deployment state, the component
registry, the master RNG and the neonatal model are caller-owned (`Model`
members). Cumulative branch diff vs `main`: ~1800+/2300− in
model+unittest — the tree is now ~600 lines smaller than before this
sweep with no behaviour change. The vector-side (string-keyed decay maps,
traps, NHH) checkpoint gaps remain as noted in §6 — gated by the
byte-identical contract.

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
              rebuilt from XML, only changeHS/changeEIR replays survive
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
