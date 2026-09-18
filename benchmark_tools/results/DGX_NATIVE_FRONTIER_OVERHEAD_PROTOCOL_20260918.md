# Prospective Native Frontier Overhead Panel

The completed21831 array demonstrated native integration, not collector
overhead. Freeze this paired experiment before any of its native outcomes.
Plan:`dgx_frontier_overhead_plan_20260918.json`, SHA-256
`58e482e6f123bc82de5e98df3c93f4bf0ca3a47d2302ffcaa47290c7d1e20685`.
The plan remains execution_authorized=false pending a verified launcher and
transferred recipe. This is not authorization of the27-run scientific panel.

## Fixed Design

Use the original four-proteome scaling input:73266 proteins,36474860 sequence
characters. Retain original input bytes/basenames/order, frozen core,
environment/runtime manifests, native arguments,20CPU/96GiB exclusive DGX
allocation, and single-task concurrency. All outputs and OrthoFinder input
copies are fresh; no search cache reuse. Input hashing and runtime checks
remain outside the worker command-wall timer and warm caches in both arms.

Eighteen sequential tasks: three methods, three adjacent paired replicates,
two observation modes. Pair block0 method order is high/satellite/OrthoFinder,
block1 satellite/OrthoFinder/high, block2 OrthoFinder/high/satellite. Arm
order is boundary/periodic in blocks0 and2 and periodic/boundary in block1.
This counterbalances order partially; with three pairs it is not perfectly
balanced. Do not change order, repetitions or dataset after outcomes.

Periodic arm uses the frozen native frontier collector. Boundary arm uses
the separate boundary control, with the same worker, pre/post observations,
completion polling and cleanup, omitting periodic reads/serialization only.
Each command retains a900-second native timeout and one-hour scheduler limit.
A timeout is a retained failure, not permission to extend that task or rerun
only the failed method. Eligibility delay is at least60seconds after transfer
and submission. Make no DGX SSH/SCP calls until all18 tasks are terminal.

## Analysis And Gates

Primary statistic: periodic worker native command wall time divided by its
paired boundary worker native command wall time, minus1. Keep all nine
pair-level values and per-method medians/ranges; no significance claim from
three pairs. Prespecified numerical budgets: each method median<=5% and
every pair<=10%. Negative differences are runtime variability, not proof
of negative observer cost. No overhead subtraction or corrected runtimes.

Validate scheduler success/no restarts, runtime/recipe/input identity,
native output partitions/pairs, complete raw observations and independent
replay. Require each native command>=60seconds. Compare canonical group
memberships and native pair sets within each method's arms, not arbitrary
group labels or output-file ordering. Differences invalidate an equivalent-
work interpretation and must remain visible. GNU-time CPU and memory fields,
cgroup memory and read durations are secondary descriptive diagnostics,
with their different accounting scopes stated.

Report numerical budget results separately from observation/environmental
validity. Preserve all original periodic interval and whole-command flags.
Boundary interval screening is unavailable, not passed. Missing or failed
pairs remain explicitly missing; do not calculate a complete-panel pass from
survivors. No flag is reclassified based on the observed overhead ratios.
Raw differences may be reported descriptively even if environmental validity
is unresolved, but not asserted to be a causal overhead bound.

## Limits

This input is the smallest real scaling workload, not the645-protein smoke,
but it does not validate eight/twelve-proteome workloads, long-run thermal
behavior, non-CPU isolation, or total instrumentation overhead. Passing a
numerical budget does not admit scientific timings. A scientific inclusion
policy and remaining isolation evidence are separate prospective gates.
