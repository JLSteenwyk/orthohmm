# Full Native Factorial Output Review

This adds a terminal output-semantics gate for the separately executed
thirteen-identity native factorial panel. It does not change the frozen
scientific pipeline, adapters, prospective run identities or resource limits.
It neither submits jobs nor authorizes the next identity on its own.

## Entry Point

`validate_native_factorial_outputs.py` requires the bound request path and
SHA-256, a fresh successful terminal scheduler observation, and an absent
output-review destination. A live or failed allocation cannot be admitted as
a completed full-native output through this entry point. The production
context is derived from the original request, plan and private-runtime
baseline, not from unverified output metadata. Live-controller and expired
controller/fresh-accounting formats are handled distinctly by the retained
terminal observer.

The expected metrics command is the wrapper's `sys.executable` plus
`sys.argv`, without the interpreter-only `-B`. Its working directory is the
frozen native core. This differs from earlier CLI-only/cached-replay commands;
those validators cannot be applied unchanged to the new wrapper.

The semantic kernel accepts an explicit context for small retained-artifact
tests. That testing interface is not production execution authorization.

## Checked Semantics

- Frozen factors and high-sensitivity search/graph/seed settings, 32 native
  CPUs, eight workers with four threads each, and exact conditional stages.
- Original and copied FASTA hashes, global unique nonempty gene IDs and
  sequences, per-species counts and the preparation ownership digest.
- Complete HMM checkpoint inventory/checksums, lexical gene order, species
  ownership, exact int32/float64 one-dimensional array shapes, valid hit
  endpoint indices and finite scores. Large arrays are memory-mapped and
  validated in bounded slices, not expanded into a dense similarity matrix.
- Complete, duplicate-free materialized/final partitions over every input
  gene, including singleton groups. Partition equality ignores group labels.
- Candidate-only and reconciliation-enabled satellite_v2 settings, candidate
  counts, seed-sidecar coverage and merge endpoints within candidate families.
- Reconciliation summary/metrics agreement, original input and tool provenance,
  zero checkpoint reuse, species-tree hash/taxon/finite-branch checks,
  complete root-HOG coverage and materialized/root partition agreement.
- Canonical, sorted, unique, cross-species native pairs with correct endpoint
  species and reported pair counts. Native pairs remain distinct from
  cluster-derived pairs and pre-clustering graph edges.
- Stable hashes before and after artifact inspection; symlinked evidence is
  rejected rather than silently followed.

The frozen pipeline lexically sorts genes before checkpoint construction.
FASTA/native species enumeration is therefore not checkpoint gene order.
Also, the frozen pipeline omits reconciliation-settings metadata entirely
when R is off; absence is explicitly checked rather than assuming an `off`
field exists. Both behaviors were confirmed against the source and all eight
retained native diagnostic cells.

For C-off/R-on, reconciliation overwrites its original clustered file. The
canonical original families are reconstructed from native root-HOG source
labels and their union of genes, then compared to the retained input-cluster
hash. For C-on/R-on, the separately retained candidate checkpoint supplies
that hash directly. Neither operation infers biological truth.

## Validation Evidence

The [retained test report](results/native_factorial_output_tests_20261004.xml)
records **436 passed**, zero failures/errors/skips, in 6.21 seconds. It covers
the new validator plus the native adapter/executor, partition and pair-format
helpers, raw resource replay/derivation, typed process/pressure review,
scheduler policy and memory-scope components.

The new tests replay all eight small native diagnostics from the committed
artifact copies using explicitly relocated temporary contexts. No HMM search,
tree inference, expensive diagnostic or production run is repeated. Mutations
test command/settings/stage/count drift, duplicated/missing/unknown genes,
re-sealed invalid checkpoint arrays/order, malformed native pairs, wrong tree
taxa, candidate seed/parameter errors, checkpoint reuse, changed evidence and
symlinks. A synthetic nonzero merge checks seed accounting and rejects a
constraint crossing candidate families; it is not native accuracy evidence.
Production-context tests exercise both actual observer result schemas with
mocked scheduler observations; they do not claim the live job has finished.

The initial test run retained four fixture-context failures (56 passed):
the test incorrectly looked for aligner metadata in R-off cells. The source
and tests were corrected to enforce the actual frozen conditional format.
Initial and corrected engineering JUnit files remain under `benchmarks/work`.

An explicit check after these additions verified all **916** helper hashes in
the repaired live plan unchanged. New source/test files are added separately;
no existing pinned benchmark Python helper is edited while job 22427 runs.
Unrelated dirty sample outputs are left untouched.

## Scope And Remaining Work

An output review reports `native_outputs_validated=true` but
`accuracy_evaluated=false`, `resource_measurements_admitted=false`,
`terminal_reviewed=false` and `next_identity_authorized=false`. A successful
scheduler state alone does not establish valid runtime, measurement,
environmental accounting or scientifically comparable results. Failed and
timed-out attempts still need their separate retained-failure/resource review.

At the latest observation during this milestone, job **22427** remained
RUNNING at 24:18 allocation elapsed (native step 23:35), with the initial HMM
search at 60.42%. No terminal outcome, cost, score or complete output review
is available. Leave this job intact. Next assemble the full terminal review
using raw resource replay, runtime before/after bindings and typed contention
evidence alongside this semantic gate; then rescore completed outputs under
the frozen benchmark protocols before any explicit next identity.

Timing measurements are shared-host observations. Competition for CPU,
memory bandwidth and I/O may affect elapsed times by an unknown and
potentially tool-dependent amount. Matching resource limits does not establish
isolated tool performance. No quiet-window/DGX prerequisite or unrelated
process/service changes are introduced. Broader independent-validation,
uncertainty, robustness, biological and publication requirements remain open.
