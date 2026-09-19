# Prospective Pressure-Enabled Inference Overhead Panel

Freeze this protocol and tested launcher before collecting native outcomes.
The new plan is [dgx_pressure_overhead_plan_20260919.json](dgx_pressure_overhead_plan_20260919.json),
SHA-256 `3950ccfa867c463ccfd7d8693dc331a3c85dc75ba67d9e8db89be4c1dda91c38`.
Its execution flag remains false until a separate recipe-bound engineering
authorization and transferred-source preflight are retained. This does not
authorize the scientific27-run panel.

## Fixed Work and Analysis

Inherit all native-work, ordering and analysis requirements from the
[original overhead protocol](DGX_NATIVE_FRONTIER_OVERHEAD_PROTOCOL_20260918.md):
four proteomes,73266proteins,36474860sequence characters; three methods;
three paired replicates per method;18sequential tasks with the same
partially counterbalanced order. Retain frozen core7f3a9e4, environment,
native arguments, input hashes/order,20CPU/96GiB exclusive allocation,
900-second native timeout and one-hour scheduler limit per task.

Primary statistic remains periodic native wall divided by paired boundary
native wall minus1. Each method median must be<=5% and every pair<=10% to
meet the numerical engineering budget; each command must be>=60seconds.
Canonical memberships/native pair sets must agree within each method's arms.
Keep all failures, missing pairs and CPU flags. No selective reruns, survivor
panel pass, post hoc threshold adjustment or overhead subtraction.

## Deliberate Changes

Both collectors now use `native_pressure=True`, validated in integration
job21868. Periodic measurements collect pressure at each observation;
boundary controls collect it only before/after the command. The worker,
completion polling, wrapper, fresh inputs and native arguments remain the
same. Pressure is diagnostic: do not infer foreign work by subtracting
host/native PSI or use its observed values to choose a new exclusion cutoff.

Every output path is relocated to
`/home/jlsteenwyk/projects/orthohmm-publication/pressure_frontier_overhead_v1`.
No old run directory, cache or failed pair is overwritten. This is a new
complete panel with changed monitoring, not selected replacement evidence
for21838. Failure-observation retention fixes are included, but frontier
identity/order failures are still failures; no topology gate is relaxed.

Use the new explicit `--plan-sha` pin when invoking
`run_dgx_frontier_overhead.py`. Authorization purpose must be
`native_pressure_frontier_incremental_overhead`, indices0..17 only,
scientific_execution_authorized=false, and bind this plan and the exact
new recipe manifest. The launcher verifies both collectors, itself, plan,
native pressure probe and parser in the recipe. Full runtime/recipe/input
hash checks remain outside native timing before and after execution.

## Scheduler Evidence and Quiet Window

Submit the array held, with concurrency1 and no requeue. Start the existing
`capture_array_scheduler.py` on the local controller host, with immutable
poll records every5seconds and a lifetime covering all18one-hour limits
plus setup allowance. Verify the collector is live and has saved its first
poll before release. Run it as a durable scheduled controller-host task,
not on DGX and not a foreground process whose lifetime ends with a turn.
If it cannot start, leave the native array held; do not begin measurements
without terminal-record collection. Collector failures must remain visible.

Complete transfers/preflight before release; require at least60seconds
between the last DGX preparation activity and native eligibility. After
release, use only local scheduler polling until all18tasks are terminal.
Do not SSH, copy files or read remote logs during that window. Capture
detailed terminal scontrol records before expiry; sacct alone is not a
replacement. Do not stop unrelated jobs or services.

## Interpretation and Verification

Numerical budget, provenance/observation validity and environmental evidence
are separate results. Boundary interval screening remains unavailable.
Unchanged CPU gates are operational screens, not quantitative causal bounds.
Pressure and frontier data are non-atomic and include wrapper activity.
Exact DGX replay is required; report any cross-runtime floating-sum
differences explicitly rather than claiming bit-identical reproduction.

This experiment estimates incremental periodic monitoring cost on the
smallest real scaling workload. It does not establish total instrumentation
cost, eight/twelve-proteome overhead, thermal stability, or absence of every
form of interference. Passing it does not automatically admit scientific
timings. Historical results remain unchanged regardless of new outcomes.

88 focused tests passed for plan derivation, all18selections, old/new plan
separation, recipe-bound pressure dependencies, collector selection, retained
native failures and controller capture. The historical launcher source is
retained as a fixture for its original manifest check. No inference task has
yet been submitted under this protocol.
