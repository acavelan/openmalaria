# Rewrite review preparation

**Status:** 32 local pull-request (PR) sets, comprising 57 focused commits.
The branches form a dependency-ordered stack starting at the original
`cfc3c8e0`. The final stack tree matches `rewrite`, including the report and
this manifest. No PRs have been opened and these review branches have not
been pushed.

Implementation reference: `70f68dbd`. Report: [REPORT.md](REPORT.md).
The detailed results and remaining limitations remain in [PERF.md](PERF.md).

## How the stack works

- `review/base` pins `cfc3c8e0`.
- Set 01 uses `review/base` as its base. **Every subsequent set uses the
  preceding row's branch as its base**, not `main`.
- Each set contains one to four commits. A set's tip is the review unit;
  intermediate commits can introduce tests before the associated implementation.
- The original implementation history is retained on `rewrite`. The new local
  series consolidates merge commits, progress-log churn and canceled experiments
  into smaller, explicitly described changes.
- Source snapshots are preserved at the selected boundaries. Historical
  research notes, census inventories and their tooling are collected in set 32
  rather than repeatedly changing throughout the implementation reviews.

This is a first review organization of the complete work. Some sets retain
intermediate representations that later sets replace—for example the
enum-selected infection models preceding the explicit model families.
The series is dependency-ordered, not a collection of independently
cherry-pickable final features.

## Review sets

The file and line counts below describe each set's diff against its immediate
base, with rename detection. They include comments, formatting, tests and data;
they are **not** the model's `cloc` code-line metric.

| Set | Local branch | Commits | Files | Added / deleted lines | Scope |
|---|---|---:|---:|---:|---|
| 01 | `review/01-safety-harness` | 3 | 60 | +16,361 / −16 | Reference fixtures, equation/RNG checks, seed override and draw-order annotations. |
| 02 | `review/02-xml-readers` | 3 | 84 | +5,500 / −5,102 | Extract scenario readers by module; enable automatic model source discovery. |
| 03 | `review/03-run-context` | 3 | 138 | +3,114 / −2,514 | Driver, clock, RNG, registries, monitoring and demographic run ownership. |
| 04 | `review/04-clinical-trees` | 2 | 10 | +795 / −809 | Decision-tree reference cases followed by the flat compiler/interpreter. |
| 05 | `review/05-formula-types` | 1 | 12 | +564 / −948 | Decay, interpolation, incidence and pathogenesis representations. |
| 06 | `review/06-clinical-drug-state` | 1 | 25 | +568 / −775 | Intermediate clinical, emergence and drug configuration/state consolidation. |
| 07 | `review/07-transmission-types` | 1 | 23 | +608 / −1,019 | Shared sampler and mosquito/transmission model machinery. |
| 08 | `review/08-infection-boundaries` | 2 | 67 | +1,684 / −2,126 | Common infection state, cold payloads and module dependency breaks. |
| 09 | `review/09-human-columns` | 3 | 41 | +887 / −340 | Population columns, infection limits, scratch storage and ordered human updates. |
| 10 | `review/10-intervention-records` | 3 | 33 | +1,358 / −1,701 | Flat schedules, component records and per-person intervention histories. |
| 11 | `review/11-repeated-work` | 2 | 18 | +138 / −86 | Reused lag vectors, empty-vaccine checks, immunity constants and related loop cleanup. |
| 12 | `review/12-xml-deduplication` | 1 | 10 | +395 / −599 | Shared reader patterns with preserved diagnostics and evaluation order. |
| 13 | `review/13-infection-equations` | 1 | 16 | +402 / −365 | Separate shared infection state from variant-specific equations. |
| 14 | `review/14-clock-file-layout` | 2 | 133 | +705 / −1,341 | Clock module, small-file folds, directory cleanup and dead accessor removal. |
| 15 | `review/15-checkpoint-module` | 2 | 126 | +1,880 / −2,612 | Dissolve the extra run-state owner; replace checkpointing and add resume cases. |
| 16 | `review/16-file-consolidation` | 2 | 149 | +2,332 / −3,423 | Wrapper/loop cleanup, related-file folds, header guards and checkpoint declarations. |
| 17 | `review/17-model-families` | 4 | 138 | +12,733 / −13,126 | Explicit infection/within-host families, common-code cleanup, isolated reflow, then drug/clinical families. |
| 18 | `review/18-xml-validation` | 1 | 6 | +116 / −139 | Shared validation, date and component-description helpers. |
| 19 | `review/19-namespace-layout` | 1 | 139 | +2,260 / −2,917 | Flatten subsystem namespaces; preserve the descriptive type names. |
| 20 | `review/20-transmission-passes` | 4 | 5 | +200 / −71 | Optimized build recipe, availability reuse, fused sweeps, batched compaction and row binding. |
| 21 | `review/21-heap-pages` | 1 | 4 | +76 / −0 | Best-effort Linux huge-page/heap policy. |
| 22 | `review/22-compaction-runs` | 1 | 1 | +39 / −10 | Move only the contiguous survivor runs after gaps. |
| 23 | `review/23-measurement-validation` | 3 | 10 | +591 / −29 | Benchmark tools, strict output comparison, reliable run capture and auditable PGO training. |
| 24 | `review/24-runtime-ownership` | 1 | 17 | +348 / −525 | Collapse the extra human-update wrapper, simplify incidence/registry ownership and intervention folds. |
| 25 | `review/25-config-decay-checkpoint` | 1 | 31 | +211 / −467 | Direct parameter storage, shared parsing, scalar decay samples, checked vector archives and EIR sums. |
| 26 | `review/26-short-if-style` | 1 | 66 | +305 / −669 | Requested compact single-statement conditionals and formatting configuration. |
| 27 | `review/27-selective-rebinding` | 1 | 2 | +6 / −9 | Rebind only population rows that actually moved. |
| 28 | `review/28-flat-species` | 1 | 15 | +114 / −57 | Contiguous host/species storage, scenario-selection resets and compaction/RNG regression coverage. |
| 29 | `review/29-hybrid-host-state` | 1 | 16 | +96 / −164 | Model-owned immunity/clinical state and its checkpoint representation. |
| 30 | `review/30-scenario-runtime` | 1 | 20 | +383 / −638 | XML/runtime initialization boundary, deferred-payload ownership and three file folds. |
| 31 | `review/31-cache-host-updates` | 1 | 8 | +47 / −42 | Drug cache hits, local non-human-host reductions, shipping schema generation and final comments. |
| 32 | `review/32-report-evidence` | 2 | 27 | Documentation/supporting data | Final report, evidence ledger, recovery guide, inventories, research notes and this manifest. |

### Large diffs that need a particular review approach

- **01: generated/reference data.** The fixtures are their own first commit.
  The RNG reference CSV alone accounts for 14,020 added lines. Inspect its
  provenance and generator separately from executable test logic.
- **02–03: cross-module boundaries.** XML extraction and run-context plumbing
  necessarily affect many callers. Each is split into three commits by the
  responsibility being moved.
- **14–16: file movement.** Use moved-code highlighting and rename detection;
  much of the diff is existing bodies acquiring a more appropriate home.
- **17: formatting is a separate commit.** Review the two semantic commits
  before it and the model-family commit after it. The reflow is not hidden
  inside the common-algorithm changes.
- **19 and 26: presentation.** Namespace qualification and short-condition
  style have broad file coverage. They are isolated from the later behavioral
  and storage changes.
- **27–31: the former large final follow-up is now five sets.** Its saved,
  measured source stages make these boundaries substantially smaller.

The short-name experiment and its reversal do not appear as separate review
commits. The rejected pyrogen fast path and its reversal are also absent from
the net series.

## Provenance and checks

### What was checked during this preparation

1. Every snapshot-based commit was constructed from its recorded source tree,
   with only the listed research/documentation paths deferred.
2. The final stack tree was compared with the complete `rewrite` tree, not just
   with `model/`. Tests, fixtures, scripts, schemas and documentation are included.
3. The flat-species, hybrid-state and XML-boundary sets were reconstructed from
   the immutable source snapshots used for the recorded measurements.
   Their tracked source bytes were checked against those snapshots.
4. `HostColumnsSuite.h`, which was still untracked when the old cumulative
   patch files were saved, was recovered from those source snapshots.
5. The pre-reflow cleanup tree in set 17 was reconstructed by a conflict-free
   merge of `488ec0ea^1` and `226b5c9c`. All 133 changed source paths compare
   identically after whitespace removal against `488ec0ea`; the subsequent
   commit applies the recorded formatting.
6. Parent relationships, commit counts and final tree identity were checked.
   The original working index was protected during reconstruction.

The complete machine-readable mapping of new commits to original source
commits/trees is stored locally at:

```text
/tmp/opencode/pr-series-20261125/series.json
/tmp/opencode/pr-series-20261125/prepare.py
```

Every generated commit also carries its source provenance in the commit
message, so the mapping does not depend solely on temporary files.

### Validation status

This preparation changes commit organization. **The regrouped commits have
not been freshly rebuilt or benchmarked.** Their final source is identical to
the implementation with the recorded plain 65/65, shipping 64/64 and strict
63/63 parity results; these remain historical results, not new test runs.

PR tips generally correspond to existing validated source boundaries.
The fixture-only and pre-reflow intermediate commits were introduced for
reviewability. Re-run the applicable build, unit/scenario, parity and resume
checks at each intended PR tip before publication; evaluate performance on
the shipping recipe where that set makes a speed claim.

The known limits still apply: code-size caps remain exceeded, the final
seven-scenario timing comparison regresses by 1.06%, and general restoration
of per-person intervention state is not covered. Regrouping changes none of
those results.

Preparation commits are unsigned: the interactive signing prompt timed out.
The repository's signing configuration was not changed.

## Review commands

For example, inspect only the flat-species set:

```sh
git log --reverse --oneline \
  review/27-selective-rebinding..review/28-flat-species

git diff --stat \
  review/27-selective-rebinding review/28-flat-species

git diff --color-moved=zebra --color-moved-ws=allow-indentation-change \
  review/27-selective-rebinding review/28-flat-species -- model unittest
```

Inspect the commits inside the model-family/formatting set:

```sh
git log --reverse --oneline \
  review/16-file-consolidation..review/17-model-families
```

Verify that the prepared series contains the complete current tree:

```sh
git diff --exit-code rewrite review/32-report-evidence
```

The stack is ready for local review and boundary refinement. Publication and
creation of the actual PRs remain a later step.
