# OpenMalaria rewrite: what changed, why, and what it achieved

**Scope:** the `rewrite` branch through implementation commit
`70f68dbd025b689f6bf3d32c3227523814a8dbb9`, compared with the original
`cfc3c8e0`. Report date: 25 November 2026.

The rewrite reorganizes the simulator, removes duplicated machinery, changes
where state is stored, and eliminates repeated work. Its central requirement
is to preserve the existing simulation results and every supported model
variant.

The measured outcome is:

- **Code size:** 18,813 → **15,598 lines of model code**, a reduction of
  **3,215 lines, or 17.1%**. Some of this is structural simplification; some
  is more compact presentation.
- **Performance:** the final executable uses **1.62%, 2.63% and 4.92% less
  user CPU time** than the optimized comparison branch at populations of
  2,000, 5,000 and 10,000. A separate seven-scenario comparison is
  **1.06% slower**. The comparison and its limits are explained below.
- **Correctness:** all **63 exact reference-output comparisons** pass in
  both the ordinary validation build and the shipping optimized build.
- **Completion:** the code-size caps and the broader performance target
  remain unmet.

This report is organized by subject. It explains the final design and the
changes that contributed to code size, performance and readability, including
trade-offs and rejected approaches. [PERF.md](PERF.md) contains the detailed
measurement ledger; [RECOVERY.md](RECOVERY.md) contains rebuilding and recovery
instructions.

Two reference names recur below: **A0** is the original executable, used for
exact-output checks; **Baxter** is the original code plus seven upstream draft
performance changes, used as the stronger speed comparison. Their precise
builds and measurement scope are given in sections 10 and 11.

## Contents

1. [The model and the rules](#1-the-model-and-the-rules)
2. [Run ownership and the XML boundary](#2-run-ownership-and-the-xml-boundary)
3. [How model variants are represented](#3-how-model-variants-are-represented)
4. [Human state and population maintenance](#4-human-state-and-population-maintenance)
5. [Interventions and shared configuration](#5-interventions-and-shared-configuration)
6. [Removing repeated simulation work](#6-removing-repeated-simulation-work)
7. [Checkpointing](#7-checkpointing)
8. [Smaller code and easier navigation](#8-smaller-code-and-easier-navigation)
9. [Compiler and operating-system optimizations](#9-compiler-and-operating-system-optimizations)
10. [Making the correctness evidence trustworthy](#10-making-the-correctness-evidence-trustworthy)
11. [Measured results](#11-measured-results)
12. [Rejected approaches and remaining work](#12-rejected-approaches-and-remaining-work)
13. [Where to read the implementation](#13-where-to-read-the-implementation)

## 1. The model and the rules

### What OpenMalaria simulates

OpenMalaria follows individual people and malaria transmission through time.
A **scenario** specifies the population, biological model choices, mosquito
transmission, treatments, interventions and requested outputs.

The main terms used throughout this report are:

| Term | Meaning here |
|---|---|
| **Host** | A simulated person who can carry malaria parasites. A **non-human host** is an animal that mosquitoes can feed on. |
| **Within-host model** | The part of the simulator describing infections, parasite density and immunity inside one person. |
| **Parasite density** | The modeled concentration of parasites, used to determine infectiousness, illness and diagnostic results. |
| **Clinical model** | The part describing illness, treatment-seeking, treatment and outcomes such as recovery or death. |
| **Pathogenesis** | How infection produces illness. |
| **Incidence** | The occurrence of new infections. The incidence calculation turns exposure into new infection counts. |
| **EIR** | *Entomological inoculation rate*: exposure to infectious mosquito bites per person over a specified time interval. |
| **Genotype** | A genetic type of parasite that the model tracks separately. |
| **Vector** | In malaria biology, a transmitting mosquito. In C++, `std::vector` means a growable array whose elements occupy consecutive memory locations. These are different meanings. |
| **Time step**, or **tick** | One simulated update interval. Depending on the model, this is commonly one or five days. |
| **Warmup** | Simulated time used to establish the population and infection state before the main study period. This differs from benchmark warmup, which means untimed preliminary executions. |
| **Calibration** | Adjustment of model initialization to the scenario's transmission target. The numerical solver used for this was excluded from this work. |
| **Monitoring** | Recording model measurements. A **survey** records them at configured dates; **continuous output** records them periodically during the run. |

The **CPU**, or central processing unit, executes the program's instructions.
The simulation is still **single-threaded**: one execution stream advances one
scenario with one random seed. Tests may run several executions, but the
production program has not become a multi-scenario or concurrent simulation
engine.

### What “the same results” requires

A **random-number generator**, abbreviated **RNG**, produces a reproducible
sequence from a **seed**, its initial input. A **draw** consumes part of that
sequence. Moving a draw earlier, skipping one, or assigning a person's random
state to somebody else changes subsequent simulation outcomes.

Similarly, **floating-point arithmetic** stores approximate real numbers with
finite precision. Mathematically equivalent expressions need not produce
identical stored values: `(a + b) + c` can differ from `a + (b + c)`.
Even replacing division by multiplication with a precomputed reciprocal can
change rounding.

Consequently, the rewrite preserves:

- Random draws, their order, their owning person and the initialization order
  that seeds them.
- The order of arithmetic operations, sums and products.
- Person traversal, equal-date intervention ordering and reporting order.
- Existing diagnostic messages and the order in which errors are detected,
  except for explicitly recorded bug fixes.
- Every infection, within-host, clinical, transmission and intervention variant.

**Bit-exact** means identical output bits, rather than merely close numerical
answers. The compiler flag `-ffp-contract=off` remains enabled: it prevents
**fused multiply-add**, a combined multiplication/addition with different
rounding from two separate operations. No fast-math relaxation was adopted.

Checkpoint file compatibility was explicitly allowed to change. Simulation
results and reproducible continuation were not relaxed. `CalcNv0`/`CalcNvO`,
the excluded mosquito-initialization solver, received no replacement, caching
or other changes in the review follow-up.

## 2. Run ownership and the XML boundary

### One owner for the state of a run

**Configuration** describes what to simulate; **state** is what changes while
the simulation runs. **Ownership** means responsibility for keeping an object
alive and releasing its resources when it is no longer needed.

Much of the original run state lived in **globals** or **statics**: storage
shared across calls rather than explicitly belonging to one simulation.
`RunContext` now groups the clock, master RNG, genotype and drug registries,
intervention schedules, monitoring state, output streams and other run-owned
data. A **registry** is a collection that resolves names or identifiers to
configured objects.

This changes the basic relationship from:

```text
model code → scattered shared state
```

to:

```text
one run → RunContext → its clock, registries, schedules and reporting
```

The extra `RunExtensions` owner was dissolved into `RunContext`. Owning raw
pointers were subsequently replaced by `std::unique_ptr`, the standard C++
single-owner pointer that automatically deletes its object. This is an
application of **RAII**, *Resource Acquisition Is Initialization*: resource
lifetime follows object lifetime, including when an error interrupts execution.

Construction and destruction remain out of line where necessary: their
definitions live in `run.cpp`, allowing `RunContext.h` to mention types without
including every implementation header.

**Contribution:** ownership is visible, manual deletion code is smaller, and
there is one main place to inspect run state. This was not initially a
code-size or speed win: the first context conversion added 269 model lines.
Passing a context also introduces pointer accesses that a global did not need.
The later simplifications reduce that machinery; they do not make every
context access free.

### Parsing is separated from execution

**XML**, *Extensible Markup Language*, is the scenario document format.
An **XSD**, *XML Schema Definition*, specifies its permitted structure.
The existing schema generator produces C++ classes representing that document;
the resulting **schema tree** is the in-memory hierarchy of parsed XML objects.
**Parsing** means turning the document text into that structured representation.

The readers were gathered into `model/xml/`, by subject: health systems,
interventions, transmission, pharmacology, monitoring and so on. Generated
schema includes are confined to that model layer, with tests allowed to use
them as well.

`initialiseScenario` now handles parsing, validation and dependency-ordered
initialization. It returns `ScenarioRuntime`, containing the population,
transmission model and the small amount of timing/reporting information the
driver needs:

```text
scenario file
    ↓
model/xml/: parse, validate, initialize
    ↓
RunContext + ScenarioRuntime
    ↓
run.cpp: warmup, calibration, simulation and output
```

The full initial schema tree is released before simulation begins.
Schema-library exception formatting also stays inside `xml/`; `main.cpp`
retains the existing exit-status behavior. Compiler dependency records confirm
that `run.cpp`, `main.cpp` and `util/CommandLine.cpp` no longer depend on generated
schema headers.

This boundary has two important details:

- A scheduled **change-EIR** action holds an `EirSeries`, ordinary time-series
  data, instead of retaining a schema object. Its deployment-time validation
  remains at deployment.
- A scheduled **health-system change** owns its copied description through an
  automatically managed action. The description is released when the action
  executes or is discarded. Conversion remains deferred because converting
  early would alter treatment/tree registration and which error occurs first.

**Contribution:** less parser knowledge in the simulation driver, explicit
deferred-object lifetime, and fewer implementation files. The XML/runtime
boundary plus three associated file folds removed 69 model lines in the final
review sequence. No standalone speed gain was established; the reason to keep
them is structure.

The initial extraction itself added 215 model lines of declarations and
conversion glue. Moving code behind a boundary was not automatically a size
reduction; reader sharing and later boundary cleanup supplied the reductions.

The separation is not complete in every type declaration: some model headers
still declare conversion functions using schema types. It is also not a
conversion of every configuration global into a plain per-run object.

### A shorter driver and fewer dependency traps

Several supporting changes make this ownership structure easier to follow:

- The simulation body moved out of `main()` into `run_simulation`.
- `Clock.h` groups the time value, conversions and simulation clock; the broad
  `Global.h` header was removed.
- Header dependency cycles were broken using small shared declarations and
  identifiers. A **dependency cycle** is an include chain in which two modules
  require each other to be defined first.
- Runtime intervention behavior moved out of an XML-only header into
  `interventions/Components.cpp`.
- An intermediate seven-stage `HostView` forwarding layer was folded back into
  one `update` function. Its steps always ran together, in one order; the
  wrapper objects did not justify a second representation of the same flow.
  Reusable EIR buffers remain explicit.
- The unit-test within-host substitute now owns its own state rather than
  consulting shared `WHMockHooks`. A **mock** is a controlled test substitute
  for a production model.

The driver/ownership/mock/intervention-fold cleanup together removed 60 model
lines. Its measured improvement was about 0.8% at both tested population sizes,
below the project's 1% minimum performance-decision threshold.

### Run initialization is more reliable, not universally isolated

The expanded `TwoRunsSuite` exposed model-selection flags and an optional
non-malaria-fever interpolator that remained set after a different scenario.
They are now reset during initialization. An **interpolator** calculates a
value between configured data points, such as an age-specific rate.

The test runs scenario1, then scenarios 9, 10, NonHumanHosts, Vivax and Penny,
then scenario1 again, and checks exact agreement. This is stronger than
repeating the same scenario twice. Mutable configuration globals still exist,
so it establishes that tested sequence, not arbitrary mixed-run isolation or
thread safety.

**Sources:** [RunContext.h](model/RunContext.h), [run.cpp](model/run.cpp),
[DocumentLoader.h](model/xml/DocumentLoader.h),
[scenario.cpp](model/xml/scenario.cpp).

## 3. How model variants are represented

### Use small formula choices and real model families differently

A **class hierarchy** represents alternatives with a base class and specialized
derived classes. A **virtual call** selects the derived implementation at
runtime. An **enum**, or enumeration, is a named set of choices; a `switch`
selects the corresponding calculation explicitly.

The rewrite tried collapsing many hierarchies into enum-selected classes.
That removed boilerplate, but a single object containing every variant's state
became too large and obscured genuinely different models. The final design
uses both representations:

| Family | Final representation and reason |
|---|---|
| Decay functions | One configured value type with named formula choices. An entire class per short formula was unnecessary. |
| Age-group interpolation | One interpolator representing the supported interpolation choices. |
| Random distribution samplers | One configured sampler with the existing distribution choices; the external numerical library still supplies its sampling algorithms. |
| Pathogenesis | Shared state and explicitly selected equations, rather than separate classes for short alternatives. |
| Mosquito emergence and transmission | Shared machinery with explicit choices. The optional simple mosquito-population-dynamics model is a member with initialization/update hooks, rather than a separate inheritance layer. |
| Infection incidence | Free functions over two values on `Human`, with model choice configured for the run. A **free function** is an ordinary function rather than a method attached to each instance. |
| Individual infections | Separate Descriptive, Dummy, Empirical, Molineaux and Penny classes, each owning its own variant-specific state. These are supported scientific model names, not interchangeable approximations introduced by the rewrite. |
| Within-host models | `WHInterface`, a shared falciparum layer, and Common, Descriptive and Vivax variants. *Falciparum* and *vivax* are different malaria parasite species. |
| Clinical models | Separate five-day decision-tree, immediate-outcome and event-scheduler classes. The event scheduler tracks explicitly timed illness and treatment events. |
| Drug models | Separate one-compartment, three-compartment and parent/metabolite conversion classes. A **compartment** is a modeled drug-distribution space; a **metabolite** is a substance produced by drug conversion. |

This avoids both extremes: many tiny classes with no independent state, and
one large class carrying the state and dispatch code for every scientific model.
**Dispatch** means choosing the implementation to execute.

The intermediate all-variant infection object grew to roughly 500 bytes.
Separating rarely used variant state brought the default descriptive object
back toward 100 bytes. The final derived classes own that state directly,
avoiding the extra allocation and pointer access of separate payload objects.
A **payload** is the variant-specific data attached to a common object.

Constructors—the functions that initialize new objects—were layered carefully
so that inheritance did not change RNG draw order. Shared bases perform no
new draws; variant initialization preserves the old draw positions.

**Contribution:** the initial twelve-family collapse saved 561 model lines
overall, but some of its representation was later replaced. Restoring the
infection/within-host class families removed another 114 lines from that
intermediate design; restoring drug/clinical classes added 46 lines before an
associated six-line parser reduction. These are recorded local changes, not
independent reductions to add to the final total.

The readability gain is a more faithful representation of the scientific
choices. It is not a claim that removing or adding virtual calls always speeds
up the program.

The four original `dynamic_cast` sites were also removed. A **downcast** asks
whether a general object is a particular specialized type. Explicit capability
interfaces, and later the relevant tagged records, replaced these runtime type
queries. Remaining checkpoint casts follow an explicit model-kind tag.

### Clinical decision trees become data plus one interpreter

A clinical **decision tree** is a sequence of questions and actions: test for
infection, select by age or probability, choose a treatment, report a result,
or deploy an intervention.

The previous implementation used 15 node classes. The replacement has:

1. A **compile step** during configuration loading, meaning conversion into an
   efficient internal representation, not generation of machine code.
2. A flat array of nodes and associated tables.
3. One **interpreter**, the function that walks those nodes for a patient.

Names are resolved to identifiers once. Probability thresholds and age
boundaries are stored in the order expected by the original model. Where
legal, consecutive identical action leaves are combined. The old
whole-tree equality/deduplication machinery was removed; identical trees may
occupy separate storage without changing execution.

This is more than replacing pointers with indices. The implementation retains
one uniform random draw at each executed random node, preserves equal-key
handling, child order, treatment/reporting behavior, and registration/error
order. Tests cover boundaries and less obvious legacy behavior.

**Contribution:** 163 model lines removed at the conversion, with one place to
read tree execution. The clinical-tree path accounted for only about 0.55% of
the representative run in that measurement, so the change was not a meaningful
whole-program speed optimization.

### Incidence and infection limits have less machinery

The per-person incidence object was reduced to `p_infected` and
`cumulative_eir_a` directly on `Human`. The equations are free functions in
`InfectionIncidenceModel.cpp`. This removes one layer of per-person allocation,
ownership and forwarding without changing the available incidence equations.

The maximum simultaneous infection count, 21, has one named owner,
`SlotAllocator`, now colocated with `WHInterface`. Here **slot** means capacity
for an infection; this is not a new memory pool.

Two existing behaviors remain distinct: importing at the limit skips creation
and its genotype draw, whereas a normal update draws a count before limiting
the accepted infections. Combining those paths would change the RNG stream.

**Sources:** [Infection.h](model/Host/WithinHost/Infection/Infection.h),
[WHInterface.h](model/Host/WithinHost/WHInterface.h),
[CMDecisionTree.h](model/Clinical/CMDecisionTree.h),
[ClinicalModel.h](model/Clinical/ClinicalModel.h),
[LSTMDrug.h](model/PkPd/LSTMDrug.h),
[InfectionIncidenceModel.cpp](model/Host/InfectionIncidenceModel.cpp).

## 4. Human state and population maintenance

### Why the final storage is hybrid

**Array of Structures**, or **AoS**, groups one object's fields together:

```text
[person 0: birth, immunity, clinical state] [person 1: ...] ...
```

**Structure of Arrays**, or **SoA**, groups the same field across people:

```text
birth:     [person 0, person 1, person 2, ...]
immunity:  [person 0, person 1, person 2, ...]
clinical:  [person 0, person 1, person 2, ...]
```

SoA can help a loop inspecting one field across many people. AoS-like grouping
can help a loop inspecting several fields of one person. The CPU's **cache**
is a small fast store of recently accessed memory; **locality** means arranging
accesses so useful nearby or recently used data can be reused.

The original human objects were more AoS-like. **The rewrite introduced
`HumanColumns`; it was not inherited from the original baseline.** Moving
everything into columns did not automatically improve every access pattern.
Many infection and clinical calculations work on one person's state at a time,
and columns containing nested vectors still lead to separately allocated data.
The initial column migration also added binding and maintenance code; its
standalone timing did not establish a speed gain.

The retained **hybrid layout** deliberately combines the two approaches:

| Data | Current owner/storage | Why |
|---|---|---|
| Birth/death information, monitoring group and deployment position | Population columns | These participate in population scans and shared bookkeeping. |
| Each person's RNG | Population column | Maintains a stable, explicitly associated stream for each person. |
| Cohort memberships and intervention/vaccine histories | Per-person containers in population columns | Shared access by monitoring and intervention code. A **cohort** is a tracked subgroup of people. |
| Per-person, per-mosquito-species parameters | One contiguous host-major table | The feeding calculation consumes those parameters together. |
| Immunity accumulators and cumulative primary-infection count | `WHInterface` object | Used together by that person's within-host calculation. |
| Last-treatment time | Clinical model | Its meaning and use belong to clinical behavior. |
| Event-scheduler clocks and previous density | `CMEventScheduler` only | Other clinical variants do not need that state. |

Moving immunity and clinical state back to their consuming objects removed
**nine population columns and two row-rebinding paths**, saving 42 model lines
in the isolated hybrid step. The model objects already move with their person,
so these fields cannot accidentally retain somebody else's column index.

### Four nested species arrays become one table

Previously each person carried four separate species vectors: raw
availability, availability, biting probability and resting probability.
**Availability** describes how accessible or attractive a host is to mosquitoes.

These now form one `HostSpecies` record per person/species pair:

```text
person 0: [species 0: four values] [species 1: four values] ...
person 1: [species 0: four values] [species 1: four values] ...
```

The entire table is one `std::vector<HostSpecies>`. **Host-major** means all
species entries for a person are adjacent before the next person's entries.
There is no separate four-vector allocation bundle for each person.

The species count is fixed for the run. A person's table segment is found from
its row number and that count. Temporary **spans**—non-owning views of
consecutive elements—are reacquired after population changes so they do not
point into old storage.

The original four sampling sites and arithmetic order are preserved.
Tests cover zero, one and three species, births, overlapping moves, gaps,
removing everyone, and exact association between each survivor and its RNG.

### Compaction moves survivors only as much as necessary

**Compaction** closes storage gaps after deaths or outmigration. **Stable
compaction** preserves survivor order. That order matters here because it
affects both person identity and subsequent calculations.

The old column-maintenance approach could repeatedly shift the tail of every
column for each removal. A first improvement collected survivor indices and
compacted each column in one batch.

However, gathering every surviving row still copied too much. Instrumentation
found approximately **2.57 billion element moves** in the then-current
10,000-person run: the operation repeated across more than 10,000 simulated
steps, even though deaths often clustered near the age-ordered population tail.

The retained algorithm groups surviving indices into consecutive runs:

```text
old rows:       A B C [removed] E F G
survivor runs:  A B C           E F G
work:          leave prefix   move this run left once
```

Simple values use `memmove`, the standard operation for copying possibly
overlapping memory ranges. Owning objects use C++ **move operations**, which
transfer their buffers rather than duplicating them. Self-moves are avoided:
moving a container onto itself can destroy the contents being preserved.
The flat species table uses the same runs, multiplied by its species count.

Finally, `Human::set_row` is called only for people whose rows actually moved.
**Rebinding** means updating a person's view of its new column row. The
unnecessary full-survivor rebinding pass was deleted.

**Contribution:** batched compaction and run-based compaction produced
historical improvements of about **4.2% and 2.1% at 10,000 people** in their
respective comparisons. Selective rebinding was below the 1% decision floor
but removes an entire redundant pass. These effects cannot simply be added:
they were measured against different intermediate trees.

### Other storage simplifications

Cohort membership changed from a per-person ordered map to a short sorted
vector of `(identifier, expiry)` pairs. A **map** associates keys with values;
the sorted vector retains the old key iteration order with less container
machinery for these small lists. Not every map was replaced: correlated
deployment samples still use per-person maps inside their column.

An ownership bug was also corrected during the restructuring: a declared
destructor had suppressed automatic move operations in the old `PerHost`
representation, so an intended move copied its nested vectors instead.
The current `PerHost` remains as a view of transmission state; it is not the
old owning container.

### What the isolated storage measurements show

These three rewrite stages used the same compiler recipe and heap policy.
Times are median user CPU seconds; measurement definitions appear in
[section 11](#11-measured-results).

| Population | Selective-rebinding baseline | Flat species | Hybrid state |
|---|---:|---:|---:|
| 2,000 | 10.57 | 10.44 | 10.31 |
| 10,000 | 44.36 | 44.01 | 42.88 |

Together, flat species and hybrid state reduced CPU time by **2.46% / 3.34%**.
Hybrid versus flat alone improved **1.25% / 2.57%**. The flat-only 10k result,
0.79%, is below the 1% decision floor.

This is direct evidence for these changes on this workload. It does not
separate reduced allocation, pointer arithmetic, locality and compiler
decisions into independent causes.

**Sources:** [HostColumns.h](model/Host/HostColumns.h),
[Human.h](model/Host/Human.h), [Population.cpp](model/Population.cpp),
[PerHost.h](model/Transmission/PerHost.h),
[layout measurements](PERF.md#isolated-layout-measurements).

## 5. Interventions and shared configuration

### Scheduled actions become records

An **intervention** changes the simulated system: a vaccine, treatment,
insecticide-treated net, indoor spraying, mosquito control, or a new health
system. A **deployment** is one scheduled application of that intervention.

The timed-deployment class tree became a flat list of `TimedDeployment`
records. Each record has a date, a named kind and a `std::variant` payload.
A C++ **variant** stores one of several permitted payload types at a time.
Human intervention components similarly use `ComponentConfig` records instead
of an eleven-class tree.

A **cursor**, `nextTimed`, identifies the next scheduled action. A final
**sentinel**, an artificial “future” record, marks the end. Equal-date actions
retain their input order through **stable sorting**: sorting that preserves the
relative order of equal keys.

Age, subgroup and availability restrictions still precede the final random
coverage decision. Moving that draw ahead of a failed restriction would consume
randomness for a person who previously consumed none.

**Contribution:** fewer allocation/dispatch/ownership paths, a visible ordered
schedule, and shared deployment code. The interventions conversion recorded
a net 221-line model reduction at that point, including removal of intervention
checkpoint state. Its sub-1% timing result was not a demonstrated speedup.

### Per-person intervention objects become history records

The per-person object families for **ITNs** (insecticide-treated nets),
**IRS** (indoor residual spraying) and **GVI** (generic vector interventions)
became `IvRecord` entries. Vaccines use `VacRecord` entries.

The entries remain in deployment-history order. Redeploying updates the
matching record without moving it. This is necessary because successive
effect multiplications must keep their original floating-point order.

Net state still distinguishes:

- Never deployed: no record.
- Active: a record with an active deployment time.
- Discarded: a retained record marked inactive.

These states cannot be merged merely because two currently have no effect:
later deployment and sampling behavior can differ. Net-hole and rip sampling,
insecticide sampling and vaccine sampling remain ordered as before.

Four duplicated intervention-effect **folds** were shared. A fold combines a
sequence into one result; here it means multiplying effects from first record
to last. The common code preserves that sequential multiplication.

Correlated deployment samples moved into population-owned storage. A
**Gaussian copula** is the existing statistical construction used to correlate
intervention probability with another quantity while retaining specified
distributions. Its per-person sample is still drawn once for a component and
reused, rather than redrawn during redeployment.

### Decay configuration is shared; only the random sample is personal

A **decay function** describes an intervention effect changing with age.
**Heterogeneity** means variation between people or deployments.

Previously, making an individual decay sample cloned an entire configured
function tree. The replacement keeps the configuration shared and stores the
sampled heterogeneity as one scalar number:

```text
before: configured decay tree → cloned tree carrying an individual sample
after:  shared decay tree + individual scalar → eval_with(age, scalar)
```

The scalar is passed through every nested node, preserving the original
formula and evaluation order. The existing sampling algorithm and draw point
remain unchanged.

**Contribution:** fewer cloned objects and less lifetime/serialization code.
It was part of the 175-line parser/checkpoint/decay cleanup; no standalone
end-to-end speed gain was established for the scalar representation.

**Sources:** [Deployments.h](model/interventions/Deployments.h),
[Components.h](model/interventions/Components.h),
[IvState.h](model/interventions/IvState.h),
[DecayFunction.h](model/util/DecayFunction.h).

## 6. Removing repeated simulation work

A **hot path** is code executed frequently enough to contribute substantially
to runtime. **Allocation churn** is repeatedly requesting and releasing heap
memory for temporary objects. The **heap** is dynamically allocated program
memory.

The retained hot-path changes remove work while preserving the sequence of
model calculations.

### Reuse temporary arrays

**Scratch storage** is temporary workspace whose contents do not represent
persistent simulation state.

- EIR arrays indexed by genotype are reused across people and time steps
  through `HostScratch`.
- `updateKappa`, which updates population infectiousness toward mosquitoes,
  reuses its per-genotype arrays and resets their values.
- The within-host infectiousness calculation reuses its imported/local
  lag-density arrays instead of allocating them for every person at every tick.

Each reused array must be overwritten or reset before its next use. An early
return can make reuse incorrect if a later consumer reads old values; this is
why the lag arrays were not casually merged with the output
`probTransGenotype_*` arrays.

A historical allocation probe found **88.2 million allocations per run**,
83.3 million of them tiny, and traced them to two lag-density vectors.
Reusing them reduced that count to **5.9 million**, about 93% fewer allocations
in that measured run. Those are historical diagnostic counts, not a fresh
allocation count for today's final executable.

### Calculate invariants once

A **loop invariant** is a value that does not change during a loop.
**Hoisting** computes it outside the repeated part.

The retained cases are:

- **Empty vaccine effects:** test once per person whether vaccine records exist.
  An empty product is exactly `1.0`, so repeated per-infection or
  per-species/per-genotype scans can be skipped.
- **Immunity survival factors:** compute host/tick constants once rather than
  repeating the same exponential calculation for each infection. Initialization
  is lazy, meaning performed only when first needed, so the old error point is
  preserved.
- **Age-dependent mosquito availability:** compute the person's age factor once
  in the transmission sweep instead of repeating it for each species.
- **Transmission row access:** acquire the person's intervention and
  availability data once per person, instead of repeatedly following the
  context/column/row path.
- **EIR totals:** return the imported, local and combined sums already calculated
  by `getEIR`, instead of summing both arrays again in the human update.

None of these licenses changing the order of additions or replacing a
mathematical operation with a differently rounded equivalent.

### Fuse the transmission sweeps

An intermediate transmission implementation filled several large temporary
arrays in one pass and reduced them in another. **Loop fusion** combines passes
when their dependencies permit it.

The final implementation uses a single host-major pass. Each accumulator
still receives contributions in ascending host order, so its arithmetic
sequence is unchanged. The intermediate staging arrays and second sweep were
removed.

**Contribution:** 29 fewer model lines and a historical 2.4% improvement
against the pre-age-hoist version, or about 1% beyond that hoist. The staged
two-pass design is not still counted as a retained optimization.

### Simplify drug caches and non-human-host updates

**PK/PD** means *pharmacokinetics/pharmacodynamics*: how drug concentration
changes in the body and how that concentration affects parasites.
The `LSTM` prefix in the drug classes identifies the Liverpool School of
Tropical Medicine model family.

The drug-response cache now performs one `find` on a hit rather than repeated
map lookups. A cache hit reuses an existing result and consumes no RNG draws.
When a new correlated parent/metabolite sample is required, its sampling order
is unchanged. This changes cache access, not the drug equations.

The mosquito update previously copied the entire non-human-host map each tick
before applying temporary reductions. It now copies one record locally,
applies that record's effects, and accumulates it. Every field's decay product
and every map-ordered sum retain their order. The decay evaluations do not
mutate model state or draw randomness, allowing this traversal change.

Together these two changes removed seven model lines. The seven-scenario
total moved by +0.16% against the preceding XML stage, below the decision floor.
They are retained for simpler access and removal of unnecessary whole-map
allocation/copying, **not for a measured standalone speed win**.

### Historical measurements of the retained mechanisms

These measurements explain why work was retained. They are not additional
current-tree benchmarks, and percentages from different rows must not be added.

| Mechanism | Recorded local result | Interpretation |
|---|---|---|
| Infectiousness lag-array reuse | 2.4% less CPU; allocations 88.2M → 5.9M | Removes a measured allocation hotspot. The early build was later identified as link-time optimized, not truly profile-trained. |
| Empty-vaccine check hoist | 1.24% less CPU | Removes repeated empty work in that workload. |
| Per-host/tick immunity constants | 1.04% less CPU | Avoids repeated expensive functions without changing infection arithmetic. |
| Transmission age-factor hoist | About 1.4–2.1%; about 1.2% in a profile-trained comparison | Scales with repeated per-person work. |
| Fused transmission pass | 2.4% against its stated earlier baseline | Removes temporary storage and a sweep. |
| Once-per-person transmission row binding | About 1.4% at 2k / 1.7% at 10k | Reduces repeated pointer/index work. |
| Returning existing EIR sums | About 0.2% at 2k / 0.1% at 10k | Kept for simpler data flow and one fewer line, below the performance floor. |

**Sources:** [Human.cpp](model/Host/Human.cpp),
[WHFalciparum.cpp](model/Host/WithinHost/WHFalciparum.cpp),
[TransmissionModel.cpp](model/Transmission/TransmissionModel.cpp),
[LSTMDrugType.cpp](model/PkPd/LSTMDrugType.cpp),
[AnophelesModel.cpp](model/Transmission/Anopheles/AnophelesModel.cpp).

## 7. Checkpointing

A **checkpoint** is a saved simulation state from which execution can resume.
**Serialization** converts that state into a file; its **wire format** is the
precise byte layout.

The old implementation spread checkpoint operators and methods through model
classes and a general-purpose stream framework. It was replaced with a
self-contained `model/checkpoint/` module:

- `Stream.h` supplies typed values and bounds-checked container operations.
- The section orchestrator defines which groups of state are saved.
- Human, transmission and remaining shared-state code live in their
  corresponding checkpoint implementation files.
- A file signature and format version reject incompatible input.
- Model headers no longer declare their own serialization chains.
- The command line uses `--save-to`, `--save-at` and `--resume`.

The current format is **version 3**. It follows the final ownership model:
the species table is contiguous, and immunity/clinical fields are saved with
their owning models. Old checkpoint-file compatibility was deliberately not
maintained.

Drug doses were also made **trivially copyable**: records with no owning
pointers or custom copy behavior. Together with an order-equivalent vector
medication queue, this allows use of the existing checked vector archive
instead of separate hand-written serialization loops.

**Contribution:** checkpoint removal and replacement saved **625 net model
lines** at that conversion, despite adding a substantial dedicated module.
The benefit comes from deleting scattered methods and generic machinery,
not from omitting the population/model/RNG state needed for continuation.

### What resume tests establish

The resume harness compares a continuous run with a run interrupted by saving
and restoring. Twelve cases cover nine intervention-free cases and three
health-system-change cases, at selected save points.

Per-person intervention deployment state remains deliberately absent from
the file. The intended checkpoint workflow is after warmup, before such
interventions have accumulated. Selected model-changing actions are replayed
on restore.

Therefore the passing resume cases establish their tested continuation
behavior. They do **not** establish general restoration after arbitrary
per-person intervention deployments.

**Sources:** [checkpoint/](model/checkpoint/), [resume.py](test/resume.py).

## 8. Smaller code and easier navigation

### Replace repeated scaffolding with data or standard operations

The **standard library** is the reusable container and algorithm collection
provided with C++. The rewrite uses it where it expresses the existing
operation directly:

| Change | Effect on code and readability |
|---|---|
| `Parameters`: two lookup maps → enum-indexed optional values | One direct representation of parameter presence/value. XML identifiers map arithmetically, including the deprecated identifier gap; missing values retain their diagnostics. |
| `ModelNameProvider` wrapper → model-name enum plus XML reader | Removes an object whose main job was forwarding a selected name. The historical header filename remains. |
| Command-line parsing | Shared option/name/print handling replaces repeated branches; the seed override is explicit and tested against the scenario seed. |
| Model-option compatibility checks | A rule table replaces repetitive checks while retaining error text and priority. |
| Monitoring callbacks | Nine similar vector-statistic callbacks become one parameterized callback; shared means and output-header construction replace duplicate code. A **callback** is a function registered for later invocation. |
| Searches and tests | `find_if` and `any_of` replace hand-written search/exists loops. |
| Deletion and maxima | `erase_if` and `max_element` express conditional removal and maximum selection directly. |
| Summation | `accumulate` replaces a local summing wrapper while preserving sequential addition. It is not an unordered parallel reduction. |
| Imported/local twin loops | Common bodies are shared while finishing one side's draws before starting the other. |
| XML readers | Local helpers share repeated date, species, availability, health-system, drug and validation handling. Error text, first-match behavior and first-error order remain pinned. |
| Small boilerplate | Unused getters, empty source files, undeclared-definition leftovers, dead fields/parameters, duplicate includes and single-caller forwarding functions were removed after caller checks. |

Some simple fields are accessed directly rather than through trivial
getters/setters. Test-only `friend` declarations—special permissions to inspect
private state—were removed where the chosen representation exposes plain
state. This trades some encapsulation for fewer wrappers; it is not a claim
that all fields should be public.

Scientific variants, live equations and validation were not deleted merely
to reduce the count. Formula-named operations and accessors that preserve a
useful storage boundary remain.

The XML reader compaction alone recorded a 238-line pass and a separate
35-line pass. Shared helpers often reduced duplication despite the absence of
large text-identical blocks; preserving error order required reading the
behavior, not just applying a duplicate-code detector.

### Merge files that describe the same concept

The file-consolidation work grouped implementations that belonged together:

- The five-day clinical chain went from six files to two.
- Three drug-variant implementation files joined `LSTMDrug.cpp`.
- Dummy infection behavior joined the common infection implementation file.
- Tiny implementations for decay, species-index validation, neonatal
  mortality, continuous output and intervention parameters joined their
  related headers where the include cost was acceptable.
- The decision-tree representation and interpreter now sit together; XML
  compilation remains in the XML layer.
- The RNG engine joined its wrapper in `random.h`, with engine expressions and
  attribution preserved.
- Treatment declarations, clock declarations, pathogenesis state and the
  infection-cap helper moved beside the concepts that use them.
- Single-file or redundant directory layers such as `core/`, `humans/`,
  `infections/`, `PkPd/Drug/` and the extra pathogenesis directory were removed
  or folded into their owners.
- The final review folded `Clinical/Episode.cpp` into `ClinicalModel.cpp`,
  and two XML-only headers into `xml/interventions.cpp`.

This was not indiscriminate header-only conversion. A **header** declares
shared C++ types/functions; an implementation file supplies bodies. Each
compiled implementation file, together with its included headers, forms a
**translation unit**. Putting a large body in a widely included header makes
many translation units read it. Such bodies stayed separate where merging
would spread dependencies or create cycles.

The model's CMake source list now discovers implementation files recursively
with change tracking, replacing a manually maintained list. **CMake** is the
tool that configures compilation and test targets. This supports file moves
without repeated build-list edits; it is not a simulation-speed change.

`Human` and `Population` consequently remain separate: one describes a person
and its state access; the other owns and maintains the collection.

**Contribution:** the principal consolidation pass reduced model files
**140 → 115** and removed **353 code lines**. Other changes followed; the
current total is **113 model files**. The navigation benefit is fewer places
to inspect one concept, not an asserted runtime gain.

### Presentation changes also count as LOC reductions

**LOC** means *lines of code*. The project's counter, `cloc`, counts physical
C++ code lines and excludes blank lines and comments. A brace on its own line
can therefore affect the number even when behavior and complexity do not change.

Three important presentation changes must be read in that light:

- `#pragma once` replaces repeated header-guard declarations in handwritten
  headers. Both prevent a header from being included more than once.
- Most subsystem **namespaces**—C++ name-grouping scopes—were flattened into
  `OpenMalaria`. This removed nested wrappers and long qualifications. Useful
  small/library namespaces remain. The associated naming-reversal/reflow wave
  recorded **488 fewer lines**; proposed shorter class names were reverted.
- A 120-column formatting pass recorded **70 fewer model lines**.
  The user-requested compact single-statement `if` style subsequently removed
  **365 model code lines** without changing expressions or control flow.
  **Reflow** means changing wrapping/spacing, not logic.

These improve consistency or reduce visual scaffolding according to the chosen
style. They do not demonstrate simpler ownership, fewer instructions or faster
execution. The 17.1% headline includes them.

**Evidence:** [PERF.md](PERF.md), rows for XML compaction, file consolidation,
namespace flattening and the single-statement `if` pass; the latter is also
described in its [dedicated section](PERF.md#single-statement-if-compaction-2026-11-25-user-requested-presentation-pass).

## 9. Compiler and operating-system optimizations

### Link-time and profile-guided optimization

An **executable build recipe** includes compiler version, flags, libraries and
training inputs, not just source code.
The measured builds use GCC, the GNU compiler collection, to compile C++.

**LTO**, *link-time optimization*, lets the compiler optimize across separately
compiled files. It can inline calls—place a function body directly at its
caller—and remove work that is only visibly redundant across file boundaries.
Release builds enable it by default.

**PGO**, *profile-guided optimization*, has two build passes:

1. Build an instrumented executable that records which paths execute.
2. Run representative scenarios, then rebuild using those observations.

The current `pgo-build.sh` trains sequentially on six scenarios:
Genotypes, Penny, VecFullTest, NonHumanHosts, MSAT_HRP2 and Vivax.
This **training corpus** is a set of workloads, rather than the single scored
example. It is the same corpus used for the final optimized comparison branch.

A major measurement error was found and corrected: an early supposed PGO
build loaded no profiles. Profile identities depend on object/generated-source
paths, and warnings about missing data had been hidden. The current recipe
reuses the build directory and treats missing/stale profiles as build errors.
Claims that PGO had stopped helping were withdrawn.

The final script also regenerates the embedded schema data after cleaning.
It builds the trained executable and schemas; unit tests are built and checked
separately in the ordinary validation build.

**Contribution:** no model LOC increase for the compiler lever. A historical
same-tree comparison measured about **1.7% from LTO** and **7.9% from LTO plus
PGO relative to the plain build**. Those numbers overlap; PGO's additional
benefit was roughly 6–7%. Exact-output checks also run on the optimized
executable, so compiler changes are not assumed to preserve bits.

### Larger operating-system memory pages

The operating system maps memory in **pages**, commonly 4 KiB.
The CPU caches virtual-to-physical address translations in a **translation
lookaside buffer**, or **TLB**. A TLB miss requires extra translation work.
**Transparent huge pages**, or **THP**, let Linux back suitable memory regions
with larger pages, commonly 2 MiB. KiB and MiB mean 1,024 and 1,048,576 bytes.

The rewrite configures heap growth through `mallopt`, the allocator's tuning
interface, and supplies `MADV_HUGEPAGE` hints for eligible writable anonymous
memory mappings. **Anonymous** means memory not backed by an ordinary file.
The hints are best-effort and do not alter stored model values.

The code increases heap growth padding, disables automatic trimming of the
unused heap tail, and repeats the hint at selected lifecycle points.
Hints on individual small column allocations were ineffective; whole eligible
heap regions were needed.

**Contribution:** this added code rather than removing it, but a historical
controlled measurement showed about **1.9% less CPU at 2k and 6.1% at 10k**.
It also changes memory-retention behavior. It is a platform optimization, not
a biological or numerical one.

This distinction matters for the final comparison: the rewrite includes these
built-in heap hints and Baxter does not. Matching environment variables is
not the same as matching heap policy. Both programs can benefit from huge-page
tuning, so their timing difference cannot all be credited to storage layout.

### Profiling changed which work was attempted

**Profiling** measures where execution spends work. Ordinary function timing
shares were misleading when LTO moved bodies into different callers. A function
name could appear expensive in one build and disappear into another in the
other build.

The investigation therefore also used:

- Allocation and library-call counters to distinguish extra calls from extra
  bookkeeping.
- Hardware **retired-instruction** samples: instructions actually completed
  by the processor.
- Temporary source counters to count population moves.
- Cache/TLB counters and **IPC**, instructions completed per processor cycle.

These exposed the compaction work and the heap-policy effect. They do not prove
that SoA alone caused a cache advantage. The old blanket “SoA wins” attribution
was corrected; the measured hybrid changes are the stronger evidence.

**Sources:** [pgo-build.sh](util/pgo-build.sh),
[CMakeLists.txt](CMakeLists.txt), [HostColumns.h](model/Host/HostColumns.h),
[main.cpp](model/main.cpp), [callcount.c](util/callcount.c).

## 10. Making the correctness evidence trustworthy

### Three kinds of checks

The project's **correctness ladder** distinguishes:

- **T1: exact-output comparisons.** Run the old and new executables on the
  same inputs and compare output bits.
- **T2: equation and behavior checks.** Compare individual equations, boundary
  cases and RNG state with fixed reference values. A **fixture** or **golden**
  is such a saved reference.
- **T3: statistical comparisons.** Compare distributions across seeds against
  calibrated tolerances if a permitted numerical change prevents exact output.
  This is prepared infrastructure, not a replacement for the passing T1 checks
  reported here.

**A0** is the original reference executable from `cfc3c8e0`, built with
`-O3 -DNDEBUG -ffp-contract=off`, GCC 16.2.1 and GSL 2.8.
GCC is the C++ compiler; **GSL** is the GNU Scientific Library, which provides
numerical and distribution routines. `-O3` enables high optimization;
`-DNDEBUG` disables debug assertions.

The reference build recipe is recorded, including a library hash. A **hash**
is a content fingerprint used to detect a different binary/library. Exact
comparison claims are tied to that recipe, rather than assumed across
arbitrary compilers and library versions.

### The harness itself was repaired

A **test harness** is the code that runs tests and decides whether they passed.
Its correctness matters as much as the simulator's.

- Binary comparison now checks actual uncompressed bytes, including row order,
  **signed zero** (`+0.0` and `-0.0`) and **NaN payloads**. NaN means
  “not a number”; different bit patterns can represent it. Merely parsing
  values into a map and comparing numbers could hide these differences.
- Continuous output is also compared where produced.
- `TwoRunsSuite` now passes the correct argument count, honors its output
  filename, and rejects missing, empty or stale output.
- Harness regressions deliberately exercise these cases.
- `--seed` was checked against the scenario's equivalent seed setting.

The source census also mapped shared state, hierarchies, XML dependencies,
random draw sites and monitoring identifiers. `// parity:` comments mark
order-sensitive code. A **parity anchor** is such a pointer to the preserved
behavior; comments do not count toward model LOC.

### Specific behavioral defects addressed

Besides the mixed-scenario initialization leakage described earlier:

- An uninitialized EIR value was fixed without changing reference outputs.
- A contradictory debug assertion in the transient-EIR path was removed;
  it aborted valid vector/change-EIR scenarios in assertion-enabled builds.
- Invalid NaN decision-tree probabilities are rejected instead of reaching
  undefined ordered-map behavior. **Undefined behavior** means C++ places no
  reliable requirements on the result; this invalid-input case was not in the
  reference scenario corpus.
- Deferred health-system lifetime and cache-hit RNG behavior now have direct
  regression coverage.

These are bounded fixes, not a claim that every inherited model issue has
been resolved.

### Final checks and their scope

**CTest** is the build system's test runner. A **box scenario** is a small
end-to-end input selected to exercise a supported model feature.

| Check | Final result | What it establishes |
|---|---:|---|
| Ordinary validation build | 65/65 CTest cases; no C++ warnings/errors | Includes unit tests, harness checks, scenarios and resume cases. |
| Shipping LTO+PGO executable | 64/64 runnable CTest cases | Harness/scenario/resume coverage for that executable; no shipping-built unit-test executable is included in this count. |
| Strict A0 parity, ordinary build | 63/63 | 53 small scenario runs plus the representative 2,000-person example with seeds 1–10. |
| Strict A0 parity, shipping build | 63/63 | Same reference comparisons with compiler optimization enabled. |
| Error-path comparisons | 6/6 | Standard output, error output and process exit status match the preserved pre-follow-up executable. |

The six error cases cover a missing file, malformed XML, absent schema,
future schema version, invalid interval and a deployment-time EIP mismatch.
**EIP**, the *extrinsic incubation period*, is the parasite's development
period inside the mosquito.

Passing these cases is substantial evidence for these inputs. It is not a
mathematical proof of correctness for all scenarios, all remaining shared
configuration, or all checkpoint/intervention combinations.

**Sources:** [parity.py](util/parity.py), [compareBin.py](util/compareBin.py),
[HarnessScriptsTest.py](unittest/HarnessScriptsTest.py),
[TwoRunsSuite.h](unittest/TwoRunsSuite.h),
[a0-recipe.json](util/a0-recipe.json).

## 11. Measured results

### Code size

The following counts use:

```sh
cloc --quiet model --include-lang="C++,C/C++ Header"
cloc --quiet model/xml --include-lang="C++,C/C++ Header"
cloc --quiet unittest --include-lang="C++,C/C++ Header"
```

Only the `code` column is the LOC measure. XML is a subset of `model`, not an
additional amount. Generated schema code, scripts and this documentation are
outside that C++ model count; test code is tracked separately.

| Scope | Reference | Current | Difference |
|---|---:|---:|---:|
| All model code | 18,813 at original baseline | **15,598 / 113 files** | **−3,215 / −17.1%** |
| XML layer | 4,889 immediately after its initial extraction | **3,558 / 10 files** | −1,331 since extraction; this is not a separate baseline saving |
| C++ unit-test code | 3,529 at original baseline | **3,683 / 30 files** | +154 |
| Final review follow-up, model | 15,718 at review entry | **15,598** | −120 and three fewer files |
| Final review follow-up, C++ tests | 3,593 at review entry | **3,683** | +90 |

The main code-reduction contributors are the checkpoint replacement,
hierarchy/dispatch simplification, XML helper sharing, file consolidation,
smaller ownership/state representations, standard-library replacements and
dead-code deletion. Namespace and formatting changes are also substantial.

There is no defensible exact percentage assigning the total saving to each
concept: passes overlap, add support code, reverse earlier representations
and merge concurrently developed changes. The local measurements throughout
this report describe their stated transformations; **the total-tree difference
is the authoritative net saving**.

| Code-size target | Current | Remaining excess |
|---|---:|---:|
| At most 14,000 model lines, including XML | 15,598 | **1,598** |
| At most 2,600 XML lines | 3,558 | **958** |

### How timing was measured

**CPU time** is processor time used by the program. The scored quantity is
**user CPU seconds**, measured by `/usr/bin/time %U`; it excludes kernel CPU
time and is not the same as elapsed, or wall-clock, time.

The final comparisons used:

- The unmodified example scenario at 2,000 people and seed 1, with separate
  population-adjusted copies for 5,000 and 10,000.
- One process at a time, pinned to CPU core 5. **Pinning** restricts where the
  operating system schedules the program.
- Alternating executable order to reduce systematic first/second-run bias.
- Two benchmark warmups followed by five measured runs per executable.
- The **median**, the middle measured value, and **MAD**, median absolute
  deviation: the median distance of samples from their median.
- The same compiler, LTO and six-scenario PGO training corpus for the final
  rewrite/Baxter comparison.

The machine was an Intel Core i7-13700K, with GCC 16.2.1 and GSL 2.8.
`GLIBC_TUNABLES`, the allocator-tuning environment variable, was unset.
The CPU frequency governor was `powersave` and
`intel_pstate/no_turbo=0` allowed turbo. These were **not** fixed-frequency or
turbo-disabled measurements, despite the stricter ideal recipe described in
older planning text.

The policy treats sub-1% micro-optimization results as insufficient grounds
to keep extra complexity; measured dispersion can require an even larger
margin. Structural simplifications may be retained below that floor.
“Below the decision floor” does not establish a statistical tie.

The benchmark script records executable hashes, scenario, core and allocator
environment, plus diagnostic memory and timing metrics. A timed-out partial
5k warmup log was preserved and excluded; complete replacement windows supply
the results below.

### Final comparison with the optimized baseline

**Baxter** is the report's label for `eval/bs-all@6603c79e`: the original
baseline plus seven upstream draft performance changes, built with LTO+PGO.
It is a stronger comparison than the original plain A0 executable.

The draft changes included temporary-array reuse, population maintenance,
availability-factor reuse, build optimization and benchmarking support.
Several were already represented in the rewrite. Their lack of additional
benefit when applied to the rewrite did not mean they were ineffective on the
original baseline; testing them there confirmed their benefit.

These are **shipped-configuration comparisons**. The compiler/training setup
matches, but the built-in huge-page/heap policy differs.

| Population | Baxter CPU seconds (MAD) | Final rewrite CPU seconds (MAD) | Rewrite CPU change |
|---|---:|---:|---:|
| 2,000 | 10.47 (.000) | **10.30 (.010)** | **−1.62%** |
| 5,000 | 23.17 (.020) | **22.56 (.010)** | **−2.63%** |
| 10,000 | 45.11 (.100) | **42.89 (.030)** | **−4.92%** |

Negative change means less CPU time. A reported MAD of `.000` reflects the
timer's 0.01-second resolution, not proven absence of variation.

**RSS**, *resident set size*, is memory currently held in physical RAM.
**RAM**, random-access memory, is the machine's working memory.
Median peak RSS for final rewrite/Baxter was 34,820/33,364 KiB at 2k,
41,068/37,136 KiB at 5k and 48,268/45,484 KiB at 10k.
Thus this comparison does not claim universally lower memory consumption.

### The broader scenario comparison still regresses

These scenarios use their original small populations. “Sum of medians” means
adding each scenario's median; it is not the median of a single combined run
and is not weighted to represent every user's workload.

| Scenario | What the label covers | Baxter CPU s (MAD) | Final CPU s (MAD) | Change |
|---|---|---:|---:|---:|
| Genotypes | Parasite genetic types | 5.47 (.000) | 5.49 (.010) | +0.37% |
| Penny | Penny infection model | 2.43 (.000) | 2.53 (.010) | +4.12% |
| MSAT_HRP2 | Mass screening and treatment using the histidine-rich-protein-2 diagnostic marker | 2.50 (.000) | 2.54 (.000) | +1.60% |
| TBV | Transmission-blocking vaccine | 2.95 (.000) | 2.99 (.000) | +1.36% |
| VecFullTest | Broad mosquito/vector-model scenario | 1.99 (.000) | 1.92 (.000) | −3.52% |
| NonHumanHosts | Mosquito feeding involving non-human hosts | 1.95 (.000) | 1.97 (.000) | +1.03% |
| Molineaux | Molineaux infection model | 1.55 (.000) | 1.60 (.000) | +3.23% |
| **Sum of scenario medians** | | **18.84** | **19.04** | **+1.06%** |

The acceptance limits were no more than 0.5% slower than Baxter at 2k/10k,
no more than 1% slower at 5k, and no more than 0.5% slower on the seven-scenario
total. **The population-size checks pass; the seven-scenario check fails.**

Penny and Molineaux are clear follow-up workloads. This window did not include
the pre-follow-up selective-rebinding executable, so it does not establish that
the new hybrid layout caused their regression.

### Original-baseline timings are historical

The last recorded same-window A0 timing ladder was:

| Population | Plain original A0 | Rewrite at that earlier snapshot | Speedup |
|---|---:|---:|---:|
| 2,000 | 12.24 s | 10.68 s | 1.146× |
| 10,000 | 55.6 s | 45.82 s | ≈1.21× |

**Speedup** is reference time divided by candidate time; it is not the same
number as percentage CPU reduction. These correspond to about 12.7% and
17.6% less CPU time.

Those are not current-tree A0 measurements. Dividing those old A0 numbers by
today's rewrite timings would mix measurement windows and produce an
unsupported headline.

**Evidence:** [final measurements](PERF.md#final-shipping-measurements),
[bench-om.sh](util/bench-om.sh). Local raw samples and machine metadata are
under `/tmp/opencode/rewrite-next-20261125/`, including
`final-benchmark-results.json` and `final-benchmark-metadata.json`.

## 12. Rejected approaches and remaining work

### Changes deliberately not retained

Rejected experiments are part of the result: they prevent the report from
crediting an attractive idea that did not survive measurement.

| Approach | Why it is absent from the final implementation |
|---|---|
| Smaller RNG wrapper, 104 → 32 bytes per row | About 0.33% benefit in its test, below the floor; it added distribution-wrapper machinery. RNG storage size alone was not the measured bottleneck. |
| Extra pyrogen-threshold and age-loop hoists | Too small to justify their complexity. A **pyrogen threshold** is a modeled density-related threshold for fever; its update includes a serial recurrence, where each iteration depends on the preceding result. |
| Pyrogen zero-numerator fast path | Slightly better in plain builds but about 1.7% slower at 10k in the shipping PGO form. The shipping build decides. |
| Run-based moves for the `humans` vector itself | No useful measured gain; the expensive work was in the columns. |
| Persistent survivor-list scratch | Slower in the tested window. This is distinct from the retained reusable compaction-run list. |
| Contiguous age-group lookup and sweep-local genotype-count hoist | Flat or slightly slower. |
| Forced/transient EIR age hoists and a no-survey monitoring guard | About ±0.1%, below the floor. |
| `-fno-plt` | Changes calls through the dynamic library linkage mechanism; mixed sub-floor results. |
| `-fipa-pta` | Enables extra compiler pointer-alias analysis, reasoning about pointers that may refer to the same data; slower in both tested sizes. |
| `-march=native` / `-mtune=native` | Selects/tunes instructions for the local CPU. An exact-output experiment was about 1.17% slower, so it was not adopted. |
| Training PGO on both 2k and 10k example runs | Helped the larger case slightly but not the primary case. The current shipped recipe uses the generic six-scenario corpus. |
| Copying GSL sampling routines into the model | Measured about 2.6% benefit but added 617 model lines. Rejected in favor of retaining the library and smaller code. |
| A broad arena or infection pool | A dedicated allocation region/pool was hypothesized, but allocation profiling found temporary lag vectors were the dominant problem. Reusing those was smaller and sufficient for that finding. |
| Forced-EIR early exit | The proposed skipped initialization work was already skipped; remaining vector work produced live monitoring values. |
| Broad vectorization or reordered reductions | **SIMD**, single instruction/multiple data, processes several values together. Existing suitable loops already vectorize; cross-person random draws and serial arithmetic constrain further exact transformations. No such numerical relaxation was used to obtain the reported result. |
| Shorter names for major classes | Reverted at the user's request. |
| More aggressive XML/exception folds | Declined where they moved validation earlier, changed diagnostics or made the shared helper longer than the original code. |
| C++ language modules | The evaluated replacement for header inclusion was not adopted; it was not needed for the retained restructuring. |

### What remains open

1. **Code-size caps:** another 1,598 model lines, including 958 XML lines,
   exceed the specified caps. Formatting alone would not resolve the structural
   part of that goal.
2. **Broader performance:** the seven-scenario total misses its allowance.
   The more ambitious original rewrite speedup aspirations are not established
   by the current evidence.
3. **Attribution:** an equal-heap-policy comparison is required before assigning
   the Baxter difference specifically to storage layout. Current measurements
   compare the complete shipped configurations.
4. **Ownership completeness:** mutable configuration globals and some
   schema-facing declarations remain.
5. **Checkpoint scope:** arbitrary post-intervention continuation is not proven
   and per-host deployment state is not serialized.
6. **Excluded work:** the mosquito-initialization solver remains outside scope;
   inherited issues and separately tracked review follow-ups are not silently
   counted as fixed.

The strongest conclusion is specific: **the model is 17.1% smaller by the
chosen physical-line metric, its ownership and representation are more
explicit, and the retained performance work improves the measured larger
example runs while preserving exact tested outputs.** It is neither complete
target achievement nor a universal speedup claim.

## 13. Where to read the implementation

| Question | Entry point |
|---|---|
| What does one run own? | [model/RunContext.h](model/RunContext.h) |
| How does a scenario become executable state? | [model/xml/scenario.cpp](model/xml/scenario.cpp), [DocumentLoader.h](model/xml/DocumentLoader.h) |
| What happens each tick? | [model/run.cpp](model/run.cpp), [Host/Human.cpp](model/Host/Human.cpp) |
| Where is human state stored and moved? | [Host/HostColumns.h](model/Host/HostColumns.h), [Population.cpp](model/Population.cpp) |
| Where are infection and immunity variants? | [Host/WithinHost/WHInterface.h](model/Host/WithinHost/WHInterface.h), [Infection/](model/Host/WithinHost/Infection/) |
| How do clinical decisions work? | [Clinical/CMDecisionTree.h](model/Clinical/CMDecisionTree.h), [xml/healthSystem.cpp](model/xml/healthSystem.cpp) |
| How are interventions scheduled and applied? | [interventions/Deployments.h](model/interventions/Deployments.h), [IvState.cpp](model/interventions/IvState.cpp), [Components.cpp](model/interventions/Components.cpp) |
| How are transmission and drug calculations organized? | [TransmissionModel.cpp](model/Transmission/TransmissionModel.cpp), [PkPd/LSTMDrug.cpp](model/PkPd/LSTMDrug.cpp) |
| What is saved and restored? | [model/checkpoint/](model/checkpoint/) |
| How are results and speed checked? | [util/parity.py](util/parity.py), [util/bench-om.sh](util/bench-om.sh), [util/pgo-build.sh](util/pgo-build.sh) |
| Where are detailed measurements and rebuild instructions? | [PERF.md](PERF.md), [RECOVERY.md](RECOVERY.md) |
