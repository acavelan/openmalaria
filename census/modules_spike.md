# C++20/23/26 Modules Spike — Decision Memo

**Question:** can C++23 modules realistically replace the header/include machinery
in the OpenMalaria rewrite (LOC KPI, modernization, simplicity) — and what does
it cost on THIS toolchain?

**Verdict: reject-for-now (bulk model/), adopt-later-partially (one narrow
slice: new phase-3 `model/core` files, optional).** The toolchain compiles
modules today, but the honest net-LOC win is **~584–795 lines (3–4%)**, the
compile-time win on this codebase is **measurably zero-to-negative**, CMake
support is **Ninja-only with two dead-end sub-features**, and everything the
modules idiom would delete is deleted anyway by the phase 3–7 redesign
(free functions, RunContext, one-module-per-domain includes). Modules are a
second build system this project does not need.

All commands run on: `gcc (GCC) 16.2.1 20260810`, `cmake 4.4.5`, Arch Linux,
worktree `/tmp/opencode/om-modules-poc` @ `51ee18ab`, branch `modules-poc`.
Spike files kept under `spike/` as evidence.

---

## 1. Baseline ("before")

* `-std=`: **C++20** (`CMakeLists.txt:34` — `set(CMAKE_CXX_STANDARD 20)`).
* Configure + Release build of `openMalaria` with gate flags
  `-O3 -DNDEBUG -ffp-contract=off`: **wall 14.41 s, user CPU 104.83 s** on 24
  cores (`/usr/bin/time` around `cmake --build … --target openMalaria -j24`).
* LOC: `cloc --quiet model --include-lang="C++,C/C++ Header" | grep ^SUM`
  → `SUM: 157 4449 9372 19054` → **19,054 code lines** (86 `.h` + 71 `.cpp`).
* `#include` directives in `model/`: **1,021** (335 in `.h`, 685 in `.cpp` —
  `spike/measure_folding.py`).
* Headers: 86 (code LOC **5,913** — the pool that would "move" to `.cppm`).
* Schema-generated headers are already confined to 13 files under `model/xml/`
  (`grep -rlE '#include *[<"](xsd/|schema/)' model | wc -l` → 13).

## 2. Toy spike — what works on gcc 16.2, exactly

**Syntax gate:** standard module syntax is *not* enabled by `-std=c++23`/`-std=c++26`:

```
$ g++ -std=c++23 -c geo.cppm
geo.cppm:2:8: note: C++20 'module' only available with '-fmodules', which is
not yet enabled with '-std=c++20'
```

The flag `-fmodules` exists (`-fmodules-ts` is its alias); it is a **non-standard
opt-in flag in GCC 16 even for C++26** (`g++ -std=c++26 -Q --help=c++` shows
`-fmodules [available in C++]`, and C++23 mode's `module` note above).

**Pattern 1 — `import std;`** (`spike/geo.cppm`, `spike/geo_impl.cpp`, `spike/main.cpp`):
works, but libstdc++ ships no prebuilt BMI on Arch — you build it once:

```
g++ -std=c++23 -fmodules-ts -c /usr/include/c++/16/bits/std.cc -o std-module.o   # 1.89 s, 31 MB std.gcm
g++ -std=c++23 -fmodules-ts -c geo.cppm    -o geo.o        # 0.126 s
g++ -std=c++23 -fmodules-ts -c geo_impl.cpp -o geo_impl.o  # 0.127 s
g++ -std=c++23 -fmodules-ts -c main.cpp    -o main.o       # 0.160 s
g++ geo.o geo_impl.o main.o -o app   # linking WITHOUT std-module.o also works
```

`-std=c++26` works identically but BMIs are **dialect-locked** — the C++23
`std.gcm` is rejected: `std: error: language dialect differs 'C++23', expected
'C++26/contracts'` → one 2.0 s rebuild per dialect. No prebuilt
`libstdc++.modules.json` exists for gcc 16 here (only a stale `15.3.0` one).

**Pattern 2 — global-fragment `#include` (the realistic bridge)**
(`spike/geo2.cppm`): compiles, links, runs. Empirical semantics confirmed:
* the GMF does **not** leak to importers (correct per standard; consumers must
  still see `std::vector` via their own include/`import std`);
* the **implementation unit does not inherit the interface's GMF** —
  `geo2_impl.cpp` needs its own `#include <vector>` or fails
  (`error: 'vector' in namespace 'std' does not name a template type`);
* outputs of the `import std` app and the GMF app are `cmp`-identical
  (`MATCH: GMF-include == import-std output`).

**Real GSL through the GMF** (`spike/gslg.cppm`, `#include <gsl/gsl_math.h>` /
`<gsl/gsl_cdf.h>` inside `module;` → `export module om.gsl;`):
`GSL GLOBAL-FRAGMENT OK`. → Blocker 4c is not a blocker.

## 3. CMake 4.4.5 integration — usable for project modules, two dead ends

Generator gate: with Unix Makefiles, generate step dies:
`Modules are supported only by Ninja, Ninja Multi-Config, and Visual Studio
generators…`. **Ninja is mandatory.** On Ninja, the working minimal pattern
(`spike/cmake/CMakeLists.txt`, verbatim, builds & runs):

```cmake
cmake_minimum_required (VERSION 3.30)
project (spike_modules CXX)
set (CMAKE_CXX_STANDARD 23)
add_compile_options (-fmodules-ts)
add_executable (app main2.cpp geo2_impl.cpp)
target_sources (app PRIVATE FILE_SET CXX_MODULES FILES geo2.cppm)
```

CMake drives gcc with its own mapper (`-fmodule-mapper=…modmap -fdeps-format=p1689r5`
scanner pass), orders interface→implementation→consumer correctly, output correct.

**Dead end A — header units.** `import <vendored/PlatformUtils.hpp>;` in a
CMake-built TU fails (`failed to read compiled module … returning to the gate`):
* gcc's dependency scanner reports **0** deps for a header-unit import
  (`g++ -Mmodules … | grep -c PlatformUtils` → `0`), so CMake never schedules
  the BMI;
* there is no `CXX_HEADER_UNIT` property in 4.4.5
  (`--help-property` → "not a CMake property"); marking the header as a plain
  `FILE_SET CXX_MODULES` member errors: *"in a FILE_SET TYPE CXX_MODULES but it
  is not scheduled for compilation."*

**Dead end B — `CXX_MODULE_STD`.** The blessed route (set
`CMAKE_CXX_STDLIB_MODULES_JSON` + target property `CXX_MODULE_STD`, gate
`CMAKE_EXPERIMENTAL_CXX_IMPORT_STD`) fails for GNU with every embedded gate
GUID found in the binary: *"requires toolchain support, but it was not
provided. Reason: Experimental `import std` support not enabled when detecting
toolchain"*. **Workaround that works:** drop a copy of `/usr/include/c++/16/bits/std.cc`
in-tree as a `FILE_SET CXX_MODULES` member (`IMPORT STD VIA CMAKE OK`, full
build 2.34 s incl. the std BMI).

**Assessment for a ~150-file model:** named project modules = mechanically OK
(Ninja + `-fmodules-ts` + file sets); but header-unit imports — the only way to
`import` the XSD/Xerces world — are unbuildable in CMake+gcc today, so every
third-party/generated header needs a hand-written wrapper `.cppm`. That is a
new artifact-per-dependency tax, not a simplification.

## 4. The real blockers, with evidence

### 4a. XSD codegen — compilable via header units, unbuildable in-tree
* The generated header is one file, `build-spike/schema/scenario.h`, **6,065
  lines**, produced by the **external `xsd cxx-tree` binary**
  (`schema/CMakeLists.txt:45-56`, `--hxx-suffix .h`). `util/SchemaTranslator`
  is a legacy Java tool **outside the build** (root `CMakeLists.txt:137` commented out).
* **Standalone gcc: header-unit import works:**
  `g++ -fmodules-ts -x c++-header -c build-spike/schema/scenario.h` → 20 MB
  `.gcm`; then `import <schema/scenario.h>;` + `scnXml::Scenario s;` compiles
  (`GENERATED-XSD HEADER-UNIT IMPORT OK`). Same for Xerces
  (`XERCES HEADER-UNIT IMPORT OK`) and both GSL spellings (`import <gsl/gsl_rng.h>;`
  and `import "gsl/gsl_rng.h";` both OK).
* **In CMake: dead** (§3 dead end A) — you'd need a custom command per header
  plus a hand-maintained mapper file outside CMake's flow.
* **Make the generator emit a module instead?** Not a template tweak — there
  are no templates (external binary, CLI flags only). The only levers are
  `--hxx-prologue*` injection or a post-gen sed that rewrites the guard to
  `export module om.scnxml;` **and** converts `scenario.cpp` to an
  implementation unit — xsdcxx's `.cpp` `#include`s its own header, which
  becomes illegal inside the module; this is a fragile codegen fork, not a tweak.
* **The honest non-blocker:** phase 2 already confined schema includes to 13
  files (`grep -rlE '#include *[<"](xsd/|schema/)' model`), all in `model/xml/`.
  A plain `#include <schema/scenario.h>` inside those files' GMFs works (§2
  pattern 2). You don't need to import the generated code at all — quarantine is
  13 files. **Viable via wrapper, dead-end as a wholesale import strategy.**

### 4b. Xerces — non-issue
**Zero** `model/` files include `xercesc` directly (`grep -rl '#include *[<"]xercesc' model | wc -l` → 0;
it enters only via xsd runtime headers). Header-unit import works standalone
(336 KB BMI); inside CMake same dead end as 4a → wrapper `.cppm` with GMF
`#include` if ever wanted.

### 4c. GSL — works, tiny quarantine
GMF pattern with real headers: `GSL GLOBAL-FRAGMENT OK` (§2). Quarantine:
**9** `model/` files include `gsl/` headers. This is the one pattern modules
handle cleanly.

### 4d. CxxTest — survives, but drags the module rig into every test binary
`unittest/` has 22 `*Suite.h` CxxTest headers driven by
`contrib/cxxtest/bin/cxxtestgen` (`unittest/CMakeLists.txt:24`). Simulation
(`spike/runner.cpp` textually `#include`s `spike/suite_shim.h`, which contains
`import geo2;`): **`TEXTUAL-INCLUDE-OF-IMPORT OK`** — an `import` declaration
inside a textually-included header is legal and compiles. So the suites don't
break *as long as* the runner TU is compiled with `-fmodules-ts` and the BMI
mapper — i.e. every unittest target inherits the module build machinery for
zero benefit. Assessment: **survivable with shim, net friction +**, and T2
gates (`ctest -R unittest`) would need re-plumbing. Keep tests include-based.

### 4e. Macro metaprogramming — nearly absent
Only **5** function-like macros defined in all of `model/`
(`grep -rhoE '^\s*#\s*define\s+\w+\(' model --include='*.h'`); the only
cross-TU one is `TRACED_EXCEPTION` (13 files) + `TRACED_EXCEPTION_DEFAULT`
(7). These live fine inside GMFs/purview. **Not a blocker.**

## 5. LOC impact — the honest number (script: `spike/measure_folding.py`)

Fold-everything ceiling (each `.h`+`.cpp` pair → one `.cppm`), code lines:

| Category | model | xml/ | Notes |
|---|---|---|---|
| 1. `.h` prototypes whose body exists in a `.cpp` (`Class::m` match) | **107** | 0 | e.g. `model/Clinical/CM5DayCommon.h` `static void init();` ← body `CM5DayCommon.cpp:50` |
| 2. include-guard scaffolding (`#ifndef`/`#define`/closing `#endif`) | **249** | 9 | 86 headers × ~3 |
| 3a. `#include` lines in `.h` (deleted outright) | **335** | 22 | |
| 3b. `#include` lines in `.cpp` | (685) | (187) | replaced ~1:1 by `import` unless domains aggregate — the aggregation win belongs to the redesign, not modules |
| 4. forward-decl lines (`class X;`) in `.h` | **104** | 1 | |
| **Deletable, modules-attributable** (1+2+3a+4) | **795 (4.2% of 19,054)** | 32 | |
| Same, counting aggressive `.cpp`-include collapsing | ≤ ~1,100 (5.8%) | | upper bound, generous |

And the accounting honesty: **5,913 code LOC live in the 86 headers**; folding
moves declarations `.h`→`.cppm` — KPI-neutral ("just moved"). The real KPI gap
is 19,054 → ≤14,000 (**−5,054 needed**). Modules-attributable deletion is at
best **~8–16% of the required delta**, and its two biggest components
(duplicate prototypes, fwd-decls) are deleted *anyway* by the phase 4b free-function
conversion (a single `export double Decay(…)` definition in a `.cppm` has no
separate prototype line — but so does a plain `double Decay(…)` in a `.cpp`
declared by nothing, which is exactly the target layout). **The LOC KPI does not
need modules.**

## 6. Perf / compile-time / FP — measured

Compile time (taskset-pinned, `TIMEFORMAT=%3R`, median of 5, warm, `-O3 -DNDEBUG -ffp-contract=off`):

| TU | include | module |
|---|---|---|
| `#include <schema/scenario.h>` + 2 std hdrs (`spike/tu_inc.cpp`) | 0.020 s | 0.005 s header-unit import (`tu_hu.cpp`) — but see §3: not buildable via CMake |
| `#include <cstdio> <vector> <cmath>` TU (`spike/plain.cpp`) | 0.145 s | **0.160 s** with `import std;` (`spike/main.cpp`) — *slower* |

Full-build reality: whole `openMalaria` is **14.4 s wall / 104.8 s user**;
average TU ≈1.5 s. The headers here are lean; on the std-heavy TU pattern that
exists in `model/`, `import std` was **measurably not faster** (31 MB gcm read
costs what preprocessing costs). There is **no compile-time case** for this
codebase.

FP/codegen: same `norm2` body compiled plain vs as module implementation unit,
`objdump -d` instruction streams **identical after normalizing symbol names**
(`INSTRUCTION STREAMS IDENTICAL (ops+FP)`; only difference is name mangling —
module-linked entities become `_ZW4geo25norm2…`). `-ffp-contract=off` keeps
applying per-TU; BMIs serialize ASTs, not machine code. So FP semantics are
safe *in principle* — **but a modules build is still a new build and the full
T1 gate (51 boxes + example@2000×10 seeds, bit-identical) applies before any
phase may claim green.**

## 7. Cost ledger on this toolchain (what "adopting" actually buys you)

1. Mandatory opt-in non-standard flag `-fmodules`/`-fmodules-ts` forever (not enabled by `-std=c++26`).
2. Mandatory Ninja generator (breaks the project's Makefile/VS default flows).
3. Dialect-locked BMIs (one std.gcm rebuild per `-std=`), 31 MB std BMI, no distro prebuilts for gcc 16.
4. Header units broken in CMake+gcc → wrapper `.cppm` per external/generated header family (xsd, xerces via xsd, gsl).
5. CMake's blessed `import std` path gated off for GNU → in-tree `std.cc` copy hack.
6. Every CxxTest runner compiled as a module consumer.
7. A second build system in the critical path of a rewrite whose gates are bit-exactness.

## 8. Verdict & minimal path

**Reject-for-now** for the existing 157 files: honest net LOC win ≈584–795
(3–4%), compile-time win ≈0 (negative on the measured pattern), and all
blockers resolve only to "wrapper shims everywhere", which is *more* machinery,
not less. Do not spend phase budget here.

**The one narrow slice that is genuinely available** (if desired at all):
*new* phase-3 `model/core` files as pure named modules — no headers exist yet,
no XSD/Xerces/CxxTest inside (GMF pattern covers any RNG/GSL needs), T1 catches
any deviation. It's cheap to *try*, not to *adopt*: if attempted, keep it
behind a separate Ninja preset, one module (`om.core`), full T1 + `ctest -R
unittest` green or revert. Skip the `import std` variant (measured slower);
use GMF includes in that module too — then modules buy literally nothing over
headers, which is the memo's point in microcosm.

**Revisit triggers:** CMake unlocks `CXX_MODULE_STD` for GNU without hacks;
gcc's scanner learns to report header-unit deps (fixes 4a/4b in-build);
`-fmodules` becomes default-on at `-std=c++20`; or a phase's redesign leaves a
demonstrable include-guard/boilerplate residue that modules alone can delete.

---
Evidence files (kept): `spike/{geo,geo2,gslg}.cppm`, `spike/geo*_impl.cpp`,
`spike/main*.cpp`, `spike/plain.cpp`, `spike/tu_{inc,hu}.cpp`,
`spike/hu_{gen,xerces,gsl,gsl2}.cpp`, `spike/{runner,suite_shim}.h/cpp`,
`spike/cmake/` (4 working/failed CMakeLists + `b/`), `spike/measure_folding.py`,
`spike/gcm.cache/` (BMIs).
