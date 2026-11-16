# OpenMalaria simplification roadmap

The objective is a smaller model implementation with explicit inputs and ownership:
ordinary data structures, equations expressed as ordinary functions, and XML confined
to initialization. Preserve numerical operation order, RNG consumption, iteration
order, output formatting, and checkpoint field order.

## Direction

The target is for `Model` to own parameters, simulation time, master RNG,
population, transmission, interventions, and reporting. Immutable constants can remain constants; mutable
scenario state must belong to a model. Do not replace globals with an active-model
pointer, thread-local state, or a service locator.

Use public structs for state. Keep methods when they clarify a local operation;
avoid accessors, wrappers, and interfaces that merely forward arguments. Dispatch
model choices at initialization or the outer update boundary where practical.
Retain polymorphism when the scientific alternatives actually benefit from it.
Do not duplicate the entire XML schema as a second configuration hierarchy.

Reporting follows the same KISS rule: pass the relevant model-owned reporting
state explicitly to clinical, transmission, or intervention code and record the
event there. Ownership does not require a separate observation pass. Keep direct
counter increments or small reporting methods where the event occurs; do not add
an event bus, observer layer, or duplicated outcome records merely to separate
reporting from calculations. Apply this consistently through all affected model
variants and callers, preserving event timing, reset behavior, and checkpoint order.

## Scope and current status

This log starts at original revision `cfc3c8e0`. The first implementation was
saved in `e8b1ba76`; subsequent milestones below describe the working tree on
`astra-refactor`, not additional commits. Keep this document updated with each
completed slice, including validation and what remains unfinished.

The application is **not yet free of global state**. In particular, the clock,
most reporting state, clinical/within-host parameters, model options, genotype
and drug tables still need ownership changes. Population age-distribution and
neonatal-mortality state are now model-owned. Tests
of independent parameter sets do not establish that two complete simulations
can safely run concurrently yet.

## Completed foundations

- `Model` owns population, master RNG, transmission, intervention definitions, and
  the parameters already migrated. It does not yet own all simulation state.
- `Human` is a small state container. PerHost is replaced by human transmission
  state and functions with explicit species parameters and intervention effects.
- Infection incidence uses state, parameters, and free functions.
- `ScenarioXml` converts incidence, human transmission, sampler, age-curve, and
  decay inputs into native C++ values.
- Age interpolation is a value-owning curve, without a hierarchy or singleton.
- Decay definitions are public data. Evaluation and lifetime sampling are free
  functions. Individuals retain only a sampled scalar; configured curves and
  samplers are no longer cloned for each individual.
- Vaccines use a vector of public `VaccineEffect` structs and free functions.
  Effects refer to model-owned definitions. Vaccine parameter registration and
  first-reported-vaccine selection no longer use globals. Efficacy evaluation
  takes an explicit date.
- Population age-distribution fitting uses local state and a native input vector.
  Maternal prevalence history and neonatal risk belong to the population. Both
  calculations take their time inputs explicitly.
- Infant mortality counters and the non-malaria parameter belong to model reporting.
  Clinical updates and survey output receive that state explicitly.
- Exact regression comparison covers uninterrupted outputs, resumed outputs,
  and decompressed checkpoint payloads against the original executable.

## Change history from the original revision

### 1. Establish explicit model ownership and remove PerHost

Saved in `e8b1ba76`.

- Added `Model`, `ModelParameters`, and a separate model-initialization entry
  point. Population and master RNG became owned values, with explicit parameters
  passed to population creation and updates.
- Reduced `Human` to personal state. Human update orchestration lives in
  `Model.cpp`; there is no separate `HumanUpdate.h` layer.
- Replaced PerHost with `Host/HumanTransmission.h/.cpp` and
  `Transmission/HostTransmissionParameters.h`. Shared mosquito-species parameters
  belong to the model; personal sampled interactions belong to the human.
  Transmission calculations receive intervention state and the effect time.
- Replaced incidence subclasses with `IncidenceKind`, `IncidenceParameters`,
  `IncidenceState`, and free functions. Moved new-infection counting into the
  model's reporting state.
- Collected personal intervention state in `HumanComponents`. Kept intervention
  updates at their original position in the per-human sequence, preserving
  random draws and the order of effects.
- Moved transmission construction implementation out of its large header into
  the matching `.cpp` file.
- Added the exact regression runner and tests for incidence isolation,
  checkpoint continuation, and human transmission behavior.

Compatibility: human transmission checkpoints retain the original four
field-major vectors despite the simpler in-memory species structs. RNG algorithms,
loop ordering, and equations were preserved rather than redesigned.

### 2. Own intervention definitions and cohort-removal rules

- Replaced static `InterventionManager` state with `Model::interventions` and
  the matching `Interventions.h/.cpp` pair.
- Made definition lookup explicit in deployments, decision trees, treatments,
  health-system initialization, monitoring setup, and checkpoint restoration.
- Moved first-event cohort-removal rules into the owned intervention state.
- Bound personal vector intervention state to its owning definition; removed
  redundant lookup/casts on redeployment.
- Placed intervention ownership before population ownership in `Model`, so the
  definitions outlive the humans referencing them.

Validation added: two intervention collections can use different identifier
maps and removal rules without sharing state. This did not yet remove the
separate vaccine registry; that is completed in stage 5.

### 3. Separate the first XML configuration boundary

- Added `ScenarioXml.h/.cpp` for XML conversion. Incidence and human-transmission
  initialization now feed ordinary native inputs into calculation functions.
- Removed generated XML types from sampler and age-curve APIs. Removed the
  accidental scenario-schema inclusion through `Global.h`/`sim.h`.
- Replaced the age-interpolation hierarchy, raw ownership, and global dummy with
  one value-owning interpolator. Maximum age is an explicit construction input.
- Preserved interpolation arithmetic, duplicate-bound behavior, validation,
  and optional sample output.
- Added tests for native configuration, interpolation boundaries, validation,
  and RNG preservation, plus a compiler dependency check without generated XML
  headers.

Validation: 53 CTest tests, including 76 unit checks; all 51 scenarios exact in
both output formats and checkpoint payloads; four additional incidence variants.
This separated selected inputs, not the entire scenario XML interface.

### 4. Replace decay-object hierarchies and individual clones

- Replaced the decay hierarchy with public `DecayFunction` data and free
  evaluation/lifetime-sampling functions. Lifetime sampling takes the step
  interval explicitly.
- Moved decay XML parsing and validation into `ScenarioXml`.
- Replaced individual decay/sampler clones with optional scalar samples across
  vaccines, nets, spraying, traps, and other vector interventions.
- Removed now-unused sampler clone APIs and simplified `SimpleDecayingValue`.
- Kept the checkpoint presence flag and sampled factor in their original order.

Validation: 53 CTest tests, including 81 unit checks; all 51 scenarios exact in
both formats. A separate comparison against the original implementation covers
all 12 decay choices, finite/non-finite evaluation, lifetime sampling, and RNG
checkpoints at one- and five-day intervals. It caught signed-NaN differences in
composed arithmetic; retaining the original arithmetic calls restored byte identity.

Reduction: 358 production lines overall, counting adapters and callers.
`DecayFunction.h/.cpp` shrank from 522 to 166 lines. Per-individual decay cloning
was eliminated; no runtime speedup has been claimed or benchmarked.

### 5. Remove vaccine globals and simplify personal vaccine state

- Deleted the global vaccine-parameter registry and first-reporting-component ID.
  The first vaccine in each model's configuration receives its reporting flag
  during initialization, preserving the original reporting rule.
- Replaced `PerHumanVaccine` and `PerEffectPerHumanVaccine` with a vector of
  `VaccineEffect` structs. `vaccinate` receives state, RNG, definition, dose limits,
  date, and genotype count. `vaccineFactor` receives state, type, date, and genotype.
- Bind each effect directly to its model-owned definition. Checkpoints still
  store the component ID, dose count, deployment date, efficacies, and sampled
  factor in the same order. Loading resolves the ID against the supplied model's
  definitions and rejects missing/non-vaccine definitions.
- Removed disabled vaccine code and obsolete global registration machinery.
- Added tests for two configurations with overlapping IDs, dose-limit RNG
  preservation, reuse of personal heterogeneity on redeployment, checkpoint
  rebinding, and explicit date/type/genotype evaluation.

Validation: all 53 CTest tests pass, including 85 unit checks. All 51 scenarios
match the original in text/binary direct and resumed outputs and checkpoint
payloads. The native vaccine header also compiles independently.

Reduction: 131 production lines overall. Vaccine XML/genotype initialization
still reads legacy configuration services; removing those is part of stage 3
of the remaining work below. Model-owned time is still outstanding: this slice
removes an adjacent registry dependency and makes vaccine calculations ready to
receive that clock explicitly.

### 6. Own demography and neonatal-mortality state

- Replaced the static `AgeStructure` class and its process-wide fitting buffers,
  parameters, cumulative proportions, and callback pointer with local fit state.
  `makeAgeStructure` accepts native age groups, maximum age, and interval, and
  returns the cumulative-proportion vector owned by `ModelParameters`.
- Passed fit context through GSL's existing `void*` callback argument. Used stack
  views for the two small input vectors and scoped ownership for the minimizer.
  Kept the optimizer, starting point, step size, stopping threshold, iteration
  limit, final parameter evaluation, and floating-point expression order.
- Moved maternal prevalence history and neonatal risk into
  `Population::neonatalMortality`. The selected diagnostic is copied into the
  model parameters, rather than retained through a process-global pointer.
- Made neonatal updates accept the population, diagnostic, step start, and
  interval explicitly. Clinical updates receive the resulting risk, retaining
  the same newborn-only random draw and `<=` comparison.
- Moved demography XML conversion and neonatal diagnostic selection/warnings to
  `ScenarioXml`. Checkpoints serialize the owned risk/history at exactly their
  original position, before the population payload.
- Added tests for concurrent versus isolated demographic fits, different time
  inputs, cumulative population rounding/order, independent five-month maternal
  histories, invalid setup inputs, and exact checkpoint continuation.

Validation: all 53 CTest tests pass, including 89 unit checks. All 51 scenarios
match the original in text/binary direct and resumed outputs and checkpoint
payloads. Native population-state headers compile without generated XML headers.
Concurrent fitting tests match isolated results; this does not establish
concurrency safety for the rest of the simulation.

Reduction: 157 production lines overall; cumulative reduction is 808 lines.
This is preparation for clock ownership, not completion of it. The outer model
still obtains dates/intervals from `sim`; the migrated calculations no longer do.
Other consumers of the diagnostic registry still use its legacy global API.

### 7. Separate the sugar-feeding equation and simplify animal-host effects

The ATSB availability solver now lives in a matching `SugarFeeding.h/.cpp` pair.
Its public interface is a numeric input struct and a free function. It has no
XML, intervention, clock, human, or model dependency. GSL solver ownership is
local and automatic. Removed the hard-coded derivative-selection macro, its
unreachable alternative, and mutable global function-pointer aliases.
The secant method, initial guess, stopping rule, iteration limit, and arithmetic
expressions are unchanged. The historical combined derivative differs from the
standalone derivative; both formulas are deliberately retained. Correcting that
is a separate model-behavior change requiring scientific review.

`advancePeriod` no longer copies the entire non-human-host map. It copies one
host's scalar state, applies its five ordered effect sequences, then accumulates
contributions in the original host order. One local reduction replaces five
near-identical loops. No intervention state, serialization, expiry rule, or RNG
sampling was changed. Oviposition multipliers remain applied sequentially;
precombining them would change floating-point association.

`AnophelesModel.cpp` shrank from 960 to 792 lines. Including both new production
files, the net reduction is 52 lines. This is the first numerical boundary in
this file, **not** completion of the transmission redesign. No speedup is claimed.

Validation includes a standalone 96-case comparison with the previous solver,
a unit test checking the resulting feeding probability, and
`test/feeding_regression.py`, which exercises overlapping sugar-bait and all five
animal-host effects against an explicitly supplied original executable. It uses
the existing larviciding fixture without duplicating that XML in the repository.
Validation passed: all 51 standard scenarios in both text and binary formats,
including direct/resumed outputs and decompressed checkpoints, match original
revision `cfc3c8e0`. The overlapping-effects scenario also matches in both formats.
The 52 non-unit CTest checks passed; the final unit run passed all 90 checks.
The first unit run used an overly strict probability-residual tolerance (1e-8);
measuring the original solver established a maximum residual of approximately
5.72e-8 in the 96-case grid. The final unit tolerance is 1e-7; the separate
characterization still compares solver output bytes exactly. Production code
was unchanged during this test correction. Final executable hashes match every
candidate manifest. Records are under `build/feeding-refactor`.

To repeat the additional overlap test:

```sh
python3 test/feeding_regression.py --build build --original /path/to/original/openMalaria --output /path/to/new-results
```

### 8. Group animal-host effects and move their XML conversion to the boundary

Replaced five parallel maps of effect vectors with one
`map<string, vector<NhhIntervention>>`. Each public record holds availability,
biting, resting-site, resting-survival, and fecundity effects. The matching
Anopheles header/implementation now provide `deploy(effect, rng, time)` and
`applyIntervention(host, effect, time)`. Application returns modified scalar host
state without mutating the deployed baseline. Both operations take time explicitly.

`ScenarioXml::readNhhIntervention` now converts and validates XML into this native
record. `TimedNonHumanHostsDeployment` only registers the converted record and
schedules/deploys it; five resize paths and the repeated XML-validation blocks
are gone. Negative availability reduction remains allowed, while killing and
fecundity reductions retain their [0,1] bounds. Existing error strings and decay
labels are preserved. SpeciesIndexChecker still rejects duplicate/missing species.
No new production files were introduced.

Effects remain in deployment-instance order; each host field retains its exact
multiplication sequence. Deployment samples the five heterogeneity factors in
the original order. Legacy checkpoints do not serialize the animal-host effects;
that behavior is deliberately unchanged, so resumed runs are compared with their
original resumed counterparts, not assumed identical to uninterrupted runs.

Added unit checks for native effect application, unchanged baseline state,
optional XML fields, bounds, and sampling order across repeated deployments.
The extra feeding regression scenario now uses heterogeneous animal-effect decay
(CV 0.4), exercising RNG consumption alongside overlapping interventions.

All 53 CTest tests passed, including 93 unit checks. All 51 standard scenarios
match original revision `cfc3c8e0` in text and binary formats, for direct outputs,
resumed outputs, and decompressed checkpoint payloads. The additional overlapping
sugar-bait/heterogeneous animal-effect scenario matches in both formats as well.
The final executable hash matches all candidate manifests. Validation artifacts
are under `build/animal-effects-refactor/final-*`.

### 9. Isolate transmission checkpoint layouts and remove duplicate field lists

The checkpoint layouts for `TransmissionModel`, `NonVectorModel`, `VectorModel`,
`AnophelesModel`, and `SimpleMPDAnophelesModel` now live in the existing
`checkpoint.cpp`. Each type has one field list used by both stream directions.
The existing virtual checkpoint entry points preserve dispatch to the selected
transmission and mosquito model. Helpers are ordinary member templates implemented
in that source file; there is no new archive framework, reflection registry,
serialization macro, version format, or production file.

Removed 103 lines from the four model headers containing inline serialization
layouts. `VectorModel.cpp` also loses its duplicate checkpoint bodies.
`checkpoint.cpp` grows from 191 to 300 lines; including that growth and all
remaining declarations, the net production reduction is 9 lines. This change
primarily improves separation of responsibilities rather than code size.

Field order, parent-before-child layout, RNG placement, omitted legacy state,
read/write dispatch and stream error handling remain unchanged. Commented-out
fields were removed rather than accidentally serialized. Model ownership and
scientific calculations were not changed.

All 53 CTest tests passed, including 93 unit checks. All 51 standard scenarios
match original revision `cfc3c8e0` byte for byte in text and binary output,
including direct/resumed execution and decompressed checkpoint payloads. The
additional stochastic overlapping-intervention scenario also matches in both
formats. A source-level comparison confirms that all five transmission field
sequences retain the original read/write order. Candidate manifest hashes match
the final executable. Records are under `build/checkpoint-layout-refactor`.

### Direction after these slices

XML conversion belongs outside model calculations because it translates and
validates an external representation. Checkpointing similarly describes a storage
format and should not obscure equations or public state. Consolidating repeated
read/write field lists is worthwhile; inventing a generic serializer solely to
hide a few field lists is not.

Moving code alone does not establish a smaller design. The last three slices
mostly clarify boundaries; substantial size reductions still require removing
representation machinery, particularly clinical decision-node subclasses and
repeated model-option selection. Measure complete production changes, including
adapters, rather than reporting only the shrinking source files.

Monitoring cannot be extracted as passive logging wholesale. For example,
`Model.cpp::summarizeHuman` can remove cohort membership after reporting, and
`Human::removeFirstEvent` flushes clinical reports and resets health-system memory.
Understand these state transitions before moving code. Prefer passing model-owned
reporting state to the existing recording sites; keep the order of counter updates
and transitions fixed. A separate observation layer is not a refactoring goal.

Interventions should remain explicit scheduled state changes plus native effect
state consumed by equations. The animal-host representation now follows that
pattern; vector-population/trap setup and daily mosquito advancement still need
the same separation. The highest architectural priority remains owned simulation
time and reporting, and the strongest candidate for larger code deletion remains
clinical decision-tree representation. Do not let small adapter extractions
become a substitute for those larger tasks.

### 10. Simplify decision-node comparisons and fix temporary-node ownership

`CMDecisionTree.cpp` shrank from 804 to 730 lines without moving its contents to
another file. Standard vector/map equality replaces repeated comparison loops;
PK/PD treatment records use C++20 defaulted equality. Branch-node comparison still
uses the original recursive equality, while child collections retain pointer
comparison. The same-object shortcut remains, including its behavior for NaNs.

`save_decision` now owns each candidate immediately. Previously, finding an
existing equivalent node returned without deleting the newly allocated candidate.
The duplicate is now destroyed. Multiple, random, and age-node builders also own
their partially constructed nodes while parsing/validation may throw. Retained
nodes still belong to the existing process-wide library; that global has **not**
been removed and the node hierarchy has **not** yet been replaced.

Preserved legacy interning rules: uncomplicated-node equality ignores memory,
cohort equality ignores component ID, and simple-treatment equality compares
blood-duration vector size but not its values. These can affect which configured
node wins and therefore execution. They are documented compatibility behavior,
not endorsed semantics. Correcting them needs a separate behavior-change decision.

Added four unit checks: duplicate-node destruction, legacy simple-treatment reuse,
legacy cohort reuse, and order-sensitive report-node equality. Running the new
checks against the pre-change implementation confirms that only the duplicate
node destruction test fails; the compatibility checks pass. After restoration
and rebuild, validation compares the final executable with original revision
`cfc3c8e0`.

All 53 CTest tests passed, including 97 unit checks. All 51 standard scenarios
match original revision `cfc3c8e0` byte for byte in text and binary formats,
including direct/resumed outputs and decompressed checkpoint payloads. Candidate
manifest hashes match the final executable. The new tests were also run against
the pre-change implementation: only duplicate destruction fails, as expected.
Records are under `build/decision-comparison-refactor`.

### Original goal and remaining plan, revisited

The objective remains a small explicit owner:

```cpp
struct Model {
    ModelParameters parameters;
    SimulationTime time;
    MasterRng masterRng;
    Population population;
    Transmission transmission;
    Interventions interventions;
    Reporting reporting;
};
```

This is an ownership sketch, not a requirement to reorder actual member lifetimes
or remove useful polymorphism. Current intervention definitions must outlive the
humans that reference them. Preserve output bytes, RNG consumption, arithmetic and
iteration order, and the historical checkpoint format throughout.

Completed: PerHost removal, the smaller Human state container, incidence/decay/
vaccine representations, several owned registries and population calculations,
initial native XML inputs, and separation of transmission checkpoint layouts.
Transmission effect bookkeeping and clinical comparison code are smaller, but the
larger subsystem redesigns are incomplete.

Outstanding: own time and reporting; own clinical/within-host configuration,
options, genotype and drug tables; replace decision-tree representation and its
global library; finish vector setup/calibration/equation/intervention boundaries;
then simplify remaining clinical, infection, within-host and per-human intervention
hierarchies case by case. Remove obsolete infrastructure as those dependencies go.

Clinical node-library ownership should be designed together with clinical root
and health-system ownership, not merely moved into a new container while global
roots continue referencing it. The comparison characterization above is groundwork
for that migration and for a compact node representation. Time/reporting ownership
remains the architectural priority; no interleaved-full-model safety claim is made.

### 11. Move vector-population and trap XML setup out of transmission classes

`ScenarioXml::initVectorIntervention` and `ScenarioXml::initVectorTrap` now perform
species lookup, validation, decay conversion, and initial parameter registration.
They populate the existing native effect/trap state directly. No parallel schema
hierarchy, extra configuration structs, or new production files were introduced.
`Interventions::init` calls these adapter functions while building its schedules.

Removed both XML-facing setup methods from `VectorModel` and both corresponding
methods from `AnophelesModel`. `AnophelesModel.cpp` no longer includes
`ScenarioXml.h`, and `VectorModel.cpp` no longer needs `SpeciesIndexChecker`.
Their runtime deployment and equation code is unchanged. Other transmission XML
dependencies, notably change-EIR behavior and broader initialization, remain.

Preserved optional-effect handling, validation order and error strings, per-species
coverage/duplicate checks, trap instance order, unnamed-trap labels, decay labels,
and the existing arrays and checkpoint layout. The legacy sugar-feeding decay
label remains `seekingDeathRateIncrease`; changing labels is outside this slice.

The complete production change removes 16 lines after including the larger XML
adapter. This is mainly a responsibility-boundary improvement: moving code is
useful when it removes input-format concerns from runtime classes and eliminates
forwarding methods, even when the total line reduction is small.

Expanded `test/feeding_regression.py` to combine seeking-death increases and
oviposition mortality with sugar bait, emergence reduction, and animal-host
interventions. The added vector effects use heterogeneous decay (CV 0.2), while
animal-host effects retain CV 0.4. Original/current comparisons cover output and
checkpoint payloads in both formats. Standard scenarios also cover vector traps.

All 53 CTest tests passed, including 97 unit checks. All 51 standard scenarios
match original revision `cfc3c8e0` byte for byte in both text and binary formats,
including direct/resumed outputs and decompressed checkpoint payloads. The expanded
stochastic overlap scenario also matches in both formats. Candidate manifest
hashes match the final executable. Records are under `build/vector-xml-refactor`.

### 12. Remove transmission.cpp and place its factories at the XML boundary

The old `Transmission/transmission.cpp` contained XML-driven construction rather
than timestep transmission behavior: mode selection, seasonality conversion,
mosquito parameters, species ordering, estimator selection, and initial EIR setup.
That code now belongs to `ScenarioXml.cpp`. Removed both `transmission.cpp` and
its forwarding `transmission.h`, plus the CMake source entry and initialization
include. `ModelInitialization` now assigns the result of `xml::createTransmission`
directly to `model.transmission`.

The public factory returns `unique_ptr<TransmissionModel>`. Internal non-vector,
vector, and mosquito factories return owning pointers as well, so exceptions
during validation do not abandon already-created model objects. Species ownership
moves into its existing container before the estimator is stored and the species
index is populated, keeping the same object addresses and construction order. Helpers are private to the XML
translation unit. Removed an empty conditional containing only commented-out code.

Preserved the numerical preparation loops, species sort/comparison, model-option
selection, diagnostic errors, RNG construction order, EIR scaling, and reporting
registration behavior. This is a responsibility/ownership change, not a model
calculation change. Non-vector construction still contains XML-dependent daily
EIR preparation, and change-EIR paths still use XML; those are explicitly unfinished.
The broad clock, monitoring, and clinical-state ownership migrations also remain.

Including the larger XML adapter, this slice removes two production files and
74 production lines. Total production code is now 31,238 lines, 1,035 fewer than
original revision `cfc3c8e0` (about 3.2%).

All 53 CTest tests passed, including 97 unit checks. All 51 standard scenarios
match original revision `cfc3c8e0` byte for byte in both text and binary formats,
including direct/resumed outputs and decompressed checkpoint payloads. The expanded
stochastic overlap scenario also matches in both formats. Candidate manifest
hashes match the final executable. Records are under `build/transmission-factory-refactor`.

### 13. Use native EIR data throughout transmission interfaces

Non-vector construction now accepts a daily `vector<double>` and an integer EIP
duration. XML conversion occurs in `ScenarioXml`. The change-EIR interface takes
those same native values; vector transmission forwards the daily values to its
species, retaining the legacy behavior of ignoring that EIP-duration argument.
The numerical preparation loops and their validation remain unchanged.

`TimedChangeEIRDeployment` retains daily values and EIP duration instead of cloning
a generated XML object and managing it with a raw pointer. Its daily buffer is
released after successful deployment, as the old XML clone was. Removed the two
XML-specific mean helpers in favor of one ordinary vector mean, preserving linear
summation order and both existing empty-input error messages at their boundaries.

There are now no direct generated `scnXml` types or schema includes in the
`model/Transmission` source/header files. This does not mean the whole include
graph is XML-free: other subsystems still have generated-code dependencies.
Transmission also still reads global time, reporting, genotype and model-option
state. Numerical initialization is not automatically a pure function merely
because it no longer accepts XML.

Added native-input unit checks for mean summation order, empty input, five-day
aggregation, the positive EIR floor, and insufficient daily values. The existing
12/12a scenarios cover non-vector change-EIR execution and checkpoint resume.
An additional forced-vector change-EIR probe aborts in the **original** executable:
`VectorModel::calculateEIR` asserts `simulationMode == forcedEIR` after the change
sets `transientEIRknown`. This is existing behavior, not a refactor regression;
it was not silently fixed. The probe and original error log are preserved under
`build/native-eir-refactor/vector-change`. Vector aggregation is tested directly
instead of claiming successful end-to-end coverage of this unsupported path.

No production files were added. Net reduction: 34 lines; cumulative reduction:
1,069 lines from the original revision.

All 53 CTest tests passed, including 99 unit checks. All 51 standard scenarios
match original revision `cfc3c8e0` byte for byte in both text and binary formats,
including direct/resumed outputs and decompressed checkpoint payloads. Candidate
manifest hashes match the final executable. The forced-vector change-EIR probe
fails in the original executable as described above; no successful end-to-end
claim is made for that path. Records are under `build/native-eir-refactor`.

### 14. Own infant-mortality reporting state

Removed the process-wide infant-death counters, intervals-at-risk counters, and
non-malaria infant-mortality parameter. `Model::reporting.infantMortality` now owns
these values in a public struct. Initialization supplies the interval count and
scenario parameter explicitly; no active-model pointer or replacement registry
was introduced. Its `allCause()` calculation uses its own vector length rather
than consulting the global timestep configuration.

Human update orchestration passes the reporting instance into clinical infant
counting. The main-simulation transition resets that instance at the existing
point. Survey output takes a const reference and threads it through both text and
binary writers, retaining the original conditional evaluation and output order.
The mortality calculation, doomed-state predicate, increment order, and
zero-exposure NaN behavior are unchanged.

Checkpoint layout remains in `checkpoint.cpp`: the owned death and risk vectors
are serialized at precisely the former global-field positions. The parameter is
still reconstructed from the scenario rather than newly serialized. Removed the
static initialization and checkpoint APIs. No production files were added;
net reduction is 21 lines (1,090 cumulatively).

Added a test interleaving two reporting states with different interval counts and
non-malaria parameters. Resetting one leaves the other unchanged. This establishes
independence of these counters/calculations, not independence of whole simulations.
Most monitoring state, callbacks, the clock, and clinical parameters remain global.

All 53 CTest tests passed, including 100 unit checks. All 51 standard scenarios
match original revision `cfc3c8e0` byte for byte in both text and binary formats,
including direct/resumed outputs and decompressed checkpoint payloads. Candidate
manifest hashes match the final executable. Records are under
`build/infant-reporting-refactor`.

### 15. Make decay XML dispatch explicit

Replaced positional string-to-enum conversion and enum-range classification with
explicit name/kind pairs and named switch cases in `ScenarioXml.cpp`. Reordering
or extending `DecayKind` no longer silently changes XML interpretation. Existing
validation, parameter calculations, recursion order, equations and checkpoint
layout are preserved. Added a test evaluating all four XML arithmetic operators.
This intentionally trades a few lines for removing an implicit ordering contract.

Clock migration was traced but not implemented in this slice: its consumers span
clinical, within-host, transmission, intervention and reporting paths. It remains
the next architectural task; no global forwarding pointer or duplicate clock was
introduced.

Validation: all 53 CTest tests passed. All 51 scenarios matched original
`cfc3c8e0` in text and binary output, direct and resumed runs, and decompressed
checkpoint payloads. Records: `build/decay-parser-review-bin` and
`build/decay-parser-review-text`. Clock ownership remains unfinished.

### 16. Remove the dummy timed deployment

Deleted `DummyTimedDeployment` and its artificial future date. Normal deployment
and checkpoint replay now check `nextTimed < timed.size()` before reading a
scheduled entry. Real deployment sorting, same-date ordering, replay selection,
RNG calls and checkpoint data are unchanged. The legacy printed `Dummy` terminator
and zero-length-scenario warning remain output formatting, without a fake object.

Added a unit test covering an empty schedule, a future deployment becoming due,
repeated calls after exhaustion, and checkpoint replay past the final deployment.
No new production files or abstractions were added. Production now totals 31,191
lines, 1,082 fewer than original revision `cfc3c8e0`; the explicit decay mapping
in the preceding slice deliberately added lines to remove hidden enum ordering.

Validation: all 53 CTest tests passed. All 51 scenarios matched original
`cfc3c8e0` in text/binary output, direct/resumed runs and decompressed checkpoint
payloads. Manifests: `build/schedule-bounds-bin` and `build/schedule-bounds-text`.
Scenario 9's `--print-interventions` stdout also matched the original exactly.

### 17. Native vector deployments and explicit animal-host XML setup

Moved animal-host deployment parsing and species validation into `ScenarioXml`.
Deployment constructors no longer modify transmission or parse animal-host XML.
Prepared animal-host parameters are stored in mosquito-species order; runtime
string-key lookups and the private duplicate parameter structure were removed.
Setup remains at each original schedule entry so validation and initialization
order are preserved, including descriptions with no deployment dates.

Consolidated four vector deployment subclasses into one small record with named
population/trap/add-animal/modify-animal actions. One explicit switch replaces
repeated dynamic casts, constructors and print dispatch. Deployment order, field
multiplication order, species traversal, error text and printed labels are retained.
Change-EIR deployments receive daily numeric data directly. Pending health-system
XML uses `unique_ptr`; clinical health-system configuration still requires XML.

Added checks for reversed species order and missing, duplicate and unknown species.
`Deployments.h` shrank from 648 to 548 lines; total production shrank by 53 lines,
including XML adapters and native data definitions, to 31,138 lines. No production
files were added. Legacy checkpoint omission of animal-host effects is unchanged.

Validation: all 53 CTest tests passed. All 51 scenarios matched original
`cfc3c8e0` text/binary direct/resumed outputs and decompressed checkpoint payloads.
The overlapping animal-host/sugar-bait regression also matched in both formats.
Records: `build/animal-xml-refactor`.

### 18. Own continuous reporting state and callback lifetimes

Moved all six continuous-reporting globals into `Model::reporting.continuous`:
registered callbacks, selected outputs, stream, stream offset, period and warmup
flag. Initialization, updates, callback registration and checkpointing receive this
state explicitly. Transmission constructors no longer register callbacks globally;
initialization registers base/vector callbacks in their former order, followed by
population callbacks. No active-model pointer or global forwarding API was added.

Output entries have named `titles` and `write` fields. The two convenience callback
forms delegate to the same registration operation. Added an isolation test using
identical output names in two registries, invoking all callback forms and destroying
one registry before continuing to use the other.

The change adds 10 production lines (31,148 total) and removes six globals, without
new production files. This does not make complete models independent yet: survey
monitoring, the clock and clinical/within-host configuration still use global state.
Reporting code continues to record events directly where they occur.

Validation: all 53 CTest tests passed. All 51 scenarios matched original
`cfc3c8e0` text/binary direct/resumed outputs (including continuous output) and
checkpoint payloads. Records: `build/continuous-state-refactor`.

### 19. Store concrete samplers as values

Log-normal and gamma parameter factories now return values. The generic factory
allocates only when choosing a sampler through the existing polymorphic interface.
Decay heterogeneity no longer allocates, copies and discards a sampler. Net hole
and rip rates, relapse timing, and required PK/PD samplers are values; optional
PK/PD rates use `optional` and retain their existing presence-based model selection.
Removed redundant special-member declarations. Sampling equations (including the
legacy variance conversion), truncation, correlated samples and RNG order remain
unchanged. Added independent-copy and concrete/polymorphic RNG-equivalence checks.

No production files were added. Production shrank by six lines to 31,142; the main
benefit is removing unnecessary ownership and allocations, not a claimed speedup.
All 53 CTest tests and all 51 exact text/binary direct/resumed/checkpoint comparisons
passed, as did the overlapping feeding regression. Records: `build/sampler-values-refactor`.

### 20. Consolidate clinical decision branches

Six branch subclasses (case type, infection origin, diagnostic, uncomplicated,
severe and cohort) now share one native representation and one execution switch.
The existing child references, factory recursion, diagnostic flag, validation and
node-sharing rules remain. Memory and cohort IDs still do not participate in
legacy equivalence; changing that would change which configured node is reused.
The other node families and global library are deliberately still visible work.

Added characterization tests before changing production code: branch kinds cannot
share nodes with each other, selected outcomes and diagnostic reporting remain,
and uncomplicated memory retains its historical sharing and validation behavior.
Duration validation respects the existing simulation-step rounding.

`CMDecisionTree.cpp` shrank from 730 to 588 lines and from 15 subclasses to 10.
No production files were added. Total production is now 31,000 lines, 1,273 fewer
than the original (about 3.9%). All 53 CTest tests passed, as did all 51 exact
text/binary direct/resumed output and checkpoint comparisons with the original.
Records: `build/clinical-branches-refactor`.

### 21. Finish the continuous-output XML boundary

Continuous reporting's stream initializer now takes its native state, an explicit
filename and the resume flag. Parsing and validating monitoring XML, and selecting
registered callbacks in XML order, live in the existing ScenarioXml adapter.
The command-line output name is supplied by main. No intermediate configuration
class or new production file was needed, and the whole-production line count is
unchanged at 31,000. Output headers, flushing, resume errors and checkpoint-offset
handling remain in the reporting implementation. A native-interface test checks
disabled output and the missing-file resume error without XML or command-line setup.

All 53 CTest tests and all 51 exact text/binary direct/resumed/checkpoint
comparisons passed. Records: `build/continuous-xml-refactor`. Reporting still reads the global simulation clock;
survey reporting ownership remains a separate outstanding migration.

### 22. Give clinical decision nodes a model lifetime

Removed the process-wide decision-node library. `ModelParameters::decisions` owns
a plain `DecisionLibrary` containing the existing vector of owning pointers.
Tree construction and health-system setup receive that library explicitly.
Intervention initialization uses the same library, and scheduled health-system
changes retain its reference so checkpoint replay and later deployments preserve
node-sharing order. The model's existing member order keeps the library alive
until after humans and intervention definitions are destroyed.

The node evaluator and equivalence rules are unchanged. A lifetime test checks
within-library sharing, cross-library independence, node destruction and use of the
surviving library after another model's parameters are destroyed. The public tree
header now forward-declares the XML input instead of including generated XML.
No production files were added; this ownership change adds 14 lines, bringing the
total to 31,014 (1,259 fewer than the original).

The active clinical root pointers, health-system parameters, diagnostics and drug
registries still use global state. This removes the global node owner, not those
remaining dependencies or the one-simulation-per-process restriction. Libraries
must outlive the trees and deployments referencing them; Model enforces that order.
All 53 CTest tests and all 51 exact text/binary direct/resumed/checkpoint
comparisons passed. Records: `build/decision-library-refactor`.

### 23. Own clinical configuration and active decision roots

The former decision library is now `ModelParameters::clinical`, a native
`ClinicalParameters` struct. It owns all four active tree roots, shared mortality
curves and episode memory, five-day care/treatment parameters, daily scheduler
parameters and clinical option flags. A single enum replaces the two model-selection
booleans. Each clinical instance refers to its model's configuration, including
humans reconstructed for checkpoint loading. Health-system changes mutate that
same configuration, so existing humans see updated values without rebinding.

Removed the corresponding clinical globals and static data members. Episode
reporting reads memory from its human's clinical model; decision-tree validation
uses the supplied configuration. Scientific parameter comments moved with the data.
Tests cover independent configurations and roots, all three model selections,
and parameter changes observed by already-created instances. Existing node-sharing
and lifetime tests still apply. No checkpoint fields were added or reordered.

Production is 30,976 lines (38 fewer in this stage), including the new native
configuration header. All 53 CTest tests and all 51 exact text/binary
direct/resumed/checkpoint comparisons passed.
Records: `build/clinical-parameters-refactor`. XML initialization is still in the
clinical implementations at this stage. Within-host treatment/diagnostic/drug
registries, model options, reporting and the simulation clock remain global.

### 24. Separate clinical XML and remove the case-management wrapper

`Clinical/ClinicalXml.h/.cpp` now provides the two clinical XML entry points:
initial configuration and health-system replacement. Native configuration is still
mutated in the original order, including decision-node creation before daily
scheduler setup. XML parsing, validation, age-curve setup and treatment registration
have moved out of the runtime model implementations. Clinical decision-node XML
factories remain in CMDecisionTree and are a separate unfinished boundary.

Removed the ESCaseManagement class/header/implementation. Its setup is folded into
the XML adapter, and its runtime operation is a small free function beside decision
tree execution. Removed its already-disabled historical test (which referred to
obsolete schema/decision-map APIs); active branch and independent-root tests remain.
No live tests were removed. The existing CMake source list now builds ClinicalXml.

Runtime implementation sizes are now 147 lines for ClinicalModel, 79 for
ImmediateOutcomes, 74 for DecisionTree5Day, and 336 for EventScheduler. The new XML
implementation is 341 lines. Including both new adapter files and deleted wrapper
files, production shrank by 56 lines to 30,920 (1,353 fewer than the original).
The ownership stage added one native data header; this extraction has no net
production-file increase. All 53 CTest tests and all 51 exact text/binary
direct/resumed/checkpoint comparisons passed.
Records: `build/clinical-xml-refactor`.

### 25. Replace the global treatment registry with values

Removed the static treatment vector, TreatmentId wrapper and add/select forwarding
APIs. Clinical configuration now owns its treatment descriptions directly. A
`Treatment` is two stage-effect durations plus triggered deployments; equivalent
regimens copy this immutable data rather than sharing a process-wide numeric ID.
Within-host code receives `const Treatment&`. The identical falciparum/vivax
methods are replaced by one nonvirtual WHInterface method that calls the selected
model's `treatSimple`, then deploys triggered interventions in the original order.

Treatment XML parsing moved to ClinicalXml and returns the native value. Duration,
stage and duplicate-effect validation and triggered-deployment construction order
are retained. Triggered deployments still borrow model-owned component definitions;
copying a treatment does not register conditions or consume random numbers again.
Tests exercise value independence, the shared dispatch path, invalid durations,
unknown stages and duplicate effects. The existing mock now uses the real shared
method and observes its virtual stage-treatment call.

Removed Treatments.cpp/.h and added Treatment.h, for one fewer production file.
Production shrank by 142 lines to 30,778 (1,495 fewer than the original). No
checkpoint fields or mathematical expressions changed. All 53 CTest tests and
all 51 exact text/binary direct/resumed/checkpoint comparisons passed. Completed
records use the `-verified` suffix under `build/treatment-values-refactor`; the
earlier runs were interrupted.
Within-host genotype, diagnostic, drug and model configuration globals remain.

### 26. Own within-host selection, immunity and body-mass configuration

`ModelParameters::withinHost` now owns the selected model kind, falciparum immunity
and maternal-protection parameters, density-history length, common-model body-mass
curve and heterogeneity, and reporting/vaccine flags. Falciparum models retain a
reference to this data, which outlives humans. Selection is one enum, replacing
six global booleans. Common infections are created and restored by explicit enum
dispatch; removed the two global factory pointers and registration wrappers from
infection initialization. The empty DummyInfection initialization API is gone.

The existing equations, sampling order, model-option priority and checkpoint field
order are retained. Moved parameter meaning comments beside the native data.
Immunity tests now configure their own data instead of modifying protected globals;
additional tests check independent existing hosts and per-model factory selection.
Infection tests construct their intended concrete type directly.

Production is 30,733 lines (45 fewer in this stage), including the new native
configuration header. All 53 CTest tests and all 51 exact text/binary
direct/resumed/checkpoint comparisons passed. Records: `build/within-host-parameters-refactor`.

This does not yet own infection-specific coefficients, latent period, genotype,
drug, diagnostic, pathogenesis or vivax configuration. Those dependencies still
prevent independent full simulations; no global forwarding context was introduced.

### 27. Separate within-host initialization from runtime equations

`WithinHostXml.h/.cpp` now owns model selection and initialization of the migrated
falciparum immunity, common-model body-mass and reporting configuration. The model
initializer calls this explicit XML entry point. Removed the corresponding static
initialization methods from WHInterface, WHFalciparum, CommonWithinHost and
DescriptiveWithinHostModel; those runtime files no longer parse XML or include
generated schema headers. Initialization bodies and their order remain unchanged.

The adapter still invokes existing infection-specific, pathogenesis, vivax and
PK/PD setup. Those subsystems retain their own globals/XML paths and are not being
claimed as migrated. Added a matching header/implementation pair; including file
headers and includes this boundary costs 40 production lines, bringing the total
to 30,773 (1,500 fewer than the original). All 53 CTest tests and all 51 exact
text/binary direct/resumed/checkpoint comparisons passed. Records: `build/within-host-xml-refactor`.

### 28. Own pathogenesis parameters and consolidate model branches

Pathogenesis parameters now belong to the model's within-host configuration.
WithinHostXml initializes them; runtime pathogenesis no longer parses XML or
uses global parameter storage. One state struct dispatches the pyrogenic,
Mueller and predetermined calculations, replacing the Submodels hierarchy.
Falciparum and vivax callers receive the owned configuration explicitly.

Native configuration tests cover all three probability branches, their historical
checkpoint field layouts, restoration, and independent non-malaria fever inputs.
Time and reporting still use existing global services; this does not establish
independent simultaneous simulations.

Commit verification: the build and all 53 CTest tests passed. All 51 scenarios
matched the original baseline in text and binary direct/resumed outputs and
checkpoint payloads. Records: `build/commit-validation-text` and
`build/commit-validation-bin`. The existing Genotypes and MSAT_HRP2 binary
direct/resumed differences remain; each mode matches its original counterpart.

## Review checkpoint: are we simplifying the right things?

What this change does: The accumulated work replaces shared registries and several
small object hierarchies with model-owned data and explicit operations. PerHost
has been removed, and native configuration is beginning to replace XML objects
inside model code.

Verdict: The direction is sound, but the original architecture goal is unfinished.
The reduction is modest (roughly 4.6% of production lines), and the larger remaining
jobs should now take priority over isolated utility cleanups. No new must-fix
regression was established in the inspected ownership/update/checkpoint paths.
The assumed execution model remains one simulation per process.

Should fix:

1. **The model is still only a partial owner** (`Model.h`, `sim.cpp`, `mon/`).
   What: The model owns population, transmission and intervention definitions.
   Problem: Current time and most reporting state still live outside it; a second
   model can change the first model's environment. Fix: migrate the clock and
   reporting state through their actual callers, with an interleaved-model test.
   If skipped: the public Model shape suggests stronger independence than exists.
   This is a known unfinished migration, not a newly introduced concurrency bug.

2. **Vector advancement mixes intervention policy with equations**
   (`Transmission/Anopheles/AnophelesModel.cpp`, `interventions/Deployments.h`).
   What: `advancePeriod` expires interventions, applies effects, solves feeding,
   and advances mosquito compartments. Problem: changing deployment bookkeeping
   requires reading numerical code; animal effects had five duplicate loops.
   Fix: the solver extraction and loop consolidation above are completed. Next,
   separate native vector inputs/XML parsing, intervention-effect evaluation, and
   daily dynamics; preserve ordered multipliers instead of aggregated products.
   If skipped: this large file remains difficult to reason about and test alone.

3. **Clinical decision trees duplicate representation machinery**
   (`Clinical/CMDecisionTree.cpp`). What: each node repeats construction, execution,
   equality checks and factory handling. Problem: the same node semantics are
   spread across those operations, with a process-wide interned-node pool.
   Fix: characterize legacy sharing/equality behavior, then replace the node
   hierarchy with compact node data and one evaluator if the total code shrinks.
   If skipped: adding a decision continues to require editing several mechanisms.
   Do not silently correct existing equality omissions during that migration.

The larger transmission task should produce three real boundaries: native setup
and calibration inputs, intervention state/effect evaluation, and equations on
mosquito state. Split by those jobs, with paired headers/implementations where
needed; simply moving member functions into more files is not the objective.
Clinical decision-tree replacement is the strongest candidate for deleting
substantial boilerplate. Clock ownership is the highest architectural priority.

Not checked: exhaustive scientific correctness, performance, cross-platform byte
identity, and full simultaneous-model execution. Existing regression agreement
cannot establish those properties.

### Cumulative production-code size

Physical lines in all `model/**/*.cpp`, `*.h`, and `*.hpp`, including comments
and new adapters; tests, generated code, documentation, and build artifacts are
excluded. Counts include every production file, not just the shrinking modules.

| Snapshot | Lines | Change from original |
| --- | ---: | ---: |
| Original `cfc3c8e0` | 32,273 | — |
| First ownership/PerHost commit `e8b1ba76` | 32,152 | −121 |
| Owned interventions and first XML boundary | 32,111 | −162 |
| Decay simplification | 31,753 | −520 |
| Vaccine ownership and state simplification | 31,622 | −651 |
| Owned demography and neonatal-mortality state | 31,465 | −808 |
| Sugar-feeding boundary and animal-effect simplification | 31,413 | −860 |
| Grouped animal effects and XML boundary | 31,411 | −862 |
| Central transmission checkpoint layouts | 31,402 | −871 |
| Decision comparisons and temporary ownership | 31,328 | −945 |
| Vector intervention/trap XML boundary | 31,312 | −961 |
| Transmission construction at XML boundary | 31,238 | −1,035 |
| Native EIR interfaces and deployment state | 31,204 | −1,069 |
| Owned infant-mortality reporting | 31,183 | −1,090 |
| Explicit decay XML dispatch | 31,206 | −1,067 |
| Bounded deployment schedule | 31,191 | −1,082 |
| Native vector deployments | 31,138 | −1,135 |
| Owned continuous reporting | 31,148 | −1,125 |
| Concrete sampler values | 31,142 | −1,131 |
| Consolidated clinical branches | 31,000 | −1,273 |
| Continuous-output XML boundary | 31,000 | −1,273 |
| Owned clinical decision nodes | 31,014 | −1,259 |
| Owned clinical configuration and roots | 30,976 | −1,297 |
| Clinical XML boundary and wrapper removal | 30,920 | −1,353 |
| Treatment values and shared dispatch | 30,778 | −1,495 |
| Within-host selection and shared parameters | 30,733 | −1,540 |
| Within-host XML boundary | 30,773 | −1,500 |

The overall reduction is still modest. The larger remaining work is concentrated
in global state, clinical decision-tree machinery, within-host configuration,
and transmission bookkeeping. Line counts alone do not measure those structural
improvements or establish a performance improvement.

## Remaining work, in implementation order

1. **Own time and make update inputs explicit.** Replace mutable `sim`/`SimData`
   state with `SimulationTime` in `Model`. Separate fixed duration conversion
   from the current step/date. Pass time through population, transmission,
   within-host, clinical, intervention, and reporting calls. Retain the existing
   beginning/end-of-step distinctions and checkpoint ordering. Test interleaved
   models with different step sizes; a global clock forwarding to `Model` is
   not completion. Demographic fitting and maternal-risk calculation now accept
   explicit time inputs; the global clock and its other consumers still remain.

2. **Finish ownership of reporting and intervention registries.** Move
   `mon::internal::runtime` under reporting. Continuous-report state and callbacks
   now belong to `Model::reporting`, with explicit registration and destruction.
   Vaccine parameter/report-component globals and infant-mortality reporting
   globals have been removed. Route reporting
   explicitly from updates. Remove process-wide registration/reset lifecycles.
   Check that sequential runs and interleaved models cannot contaminate each other.

3. **Finish within-host configuration and remaining initialization dependencies.**
   Clinical parameters, roots, treatment values and shared falciparum/body-mass
   configuration now belong to the model. Clinical and shared within-host setup
   have explicit XML adapters. Move remaining infection-specific and vivax parameters,
   option flags, treatment tables, genotype/drug data, and age curves into native
   configuration owned by the model. Parse XML once at the boundary. Keep mutable
   health-system changes explicit. Remove parallel option booleans and factories
   that rediscover the selected model for each human or infection.

4. **Simplify clinical decision trees.** `CMDecisionTree.cpp` repeats constructors,
   execution, equality/dynamic-cast checks, and factory code across many node types.
   Explore a compact node vector with child indices and one evaluator. Preserve
   action ordering, random draws, shared-node behavior, and validation. Measure
   the complete replacement, including XML conversion, before adopting it.
   Repeated equality loops have now been removed and duplicate-node leaks fixed.
   Six branch subclasses now share one representation and evaluator.
   Characterization tests preserve legacy ignored cohort IDs, memory and
   blood-duration values. The node library now belongs to model parameters.
   Active clinical roots and health-system configuration now belong to the model.
   Ten remaining subclasses and the decision-tree XML factories are candidates
   for further representation and initialization work.

5. **Reduce the remaining per-human hierarchies case by case.** Clinical models,
   within-host models, infections, and vector intervention state still use runtime
   polymorphism. Separate shared parameters from individual state first. Then
   choose the smallest representation for each family: a struct and enum, a
   variant, or retained polymorphism. Preserve RNG identity and ordering when
   changing containers; do not batch/reorder humans as an incidental optimization.

6. **Simplify transmission initialization and repeated vector bookkeeping.**
   `AnophelesModel.cpp` and `CalcNvO.cpp` remain substantial. Separate XML setup,
   calibration, and daily equations. Consolidate repeated intervention-effect
   bookkeeping where the operations really are identical. Preserve floating-point
   reduction order and historical rounding. Vector-population and trap XML
   setup now lives in ScenarioXml; calibration, remaining XML paths, and the
   full separation of daily dynamics from intervention evaluation remain. The
   old transmission.cpp/.h factory pair has been removed; construction is now
   exposed through ScenarioXml with owning pointer results. Non-vector and
   change-EIR APIs now take native daily values; direct XML types have been
   removed from transmission sources and headers.

7. **Remove obsolete infrastructure as dependencies disappear.** Delete unused
   factories, clone methods, global initialization/reset APIs, raw owning pointers,
   forwarding headers, disabled code, and redundant checkpoint adapters. Make
   includes explicit and keep header/implementation pairs coherent. Avoid creating
   a generic framework to replace a handful of equations.

## Acceptance for each slice

- Count the entire production-code delta, including adapters and callers.
- Test the affected equations, validation, and stochastic branches independently.
- Run the unit and scenario tests, then exact comparisons on the same toolchain.
- Compare each execution mode with its original counterpart: `Genotypes` and
  `MSAT_HRP2` already differ between direct and resumed binary output.
- Benchmark allocation-heavy or frequently executed paths when claiming speedups.
- State remaining globals/XML dependencies honestly. Cross-platform byte identity
  is a separate question; same-toolchain regression tests do not establish it.


## Reproducing the validation

Build the original revision and candidate with the same compiler, dependencies,
and CMake flags. Keep an isolated original executable; do not rebuild the baseline
from the refactored source or replace expected fixtures to make a comparison pass.
The checked-in `test/exact_regression.py` can capture and compare both executables:

```sh
python3 test/exact_regression.py --build build --executable /path/to/original/openMalaria --output /path/to/baseline-bin --format bin --jobs 4
python3 test/exact_regression.py --build build --output /path/to/candidate-bin --baseline /path/to/baseline-bin --format bin --jobs 4
```

Repeat with `--format text` and separate output directories, and run
`ctest --test-dir build --output-on-failure`. The runner requires fresh output
directories and records hashes in `manifest.json`. It compares uncompressed output
bytes and decompressed checkpoint payloads, excluding gzip container metadata.

Local validation records are under `build/refactor-validation`,
`build/decay-refactor`, `build/vaccine-refactor`, `build/population-refactor`, `build/feeding-refactor`, `build/animal-effects-refactor`, `build/checkpoint-layout-refactor`,
`build/decision-comparison-refactor`, `build/vector-xml-refactor`, and
`build/transmission-factory-refactor`, `build/native-eir-refactor`, and
`build/infant-reporting-refactor`. These build artifacts are
not committed; the stage descriptions and acceptance criteria in this document
are the durable record. Record final verification for every subsequent stage.
