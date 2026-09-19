# Prospective Lineage Collector Overhead Panel

## Scope

Freeze before collecting new paired timings. Native diagnostics 21995-21997
passed provenance, raw replay and output equivalence, but retained 1/7/0
narrow CPU flags. These outcomes justify testing the functioning collector's
cost, not declaring the DGX free of interference. Preserve the earlier
failed panels, all flags and the limitations in `LINEAGE_NATIVE_RESULT_21995.md`.

## Design

Run a complete fresh 18-task panel: three methods, three paired repetitions,
two collector modes. Derive exact commands, order, inputs, seeds, timeout,
resource limits and output semantics from pressure-overhead-v2 plan SHA-256
`b644e165dbf4d0beabf1cf4d9b6c314de522e3ebd1b91598ebebea99094c8fff`.
Relocate only output/cache paths to `lineage_collector_overhead_v1`.
The workloads remain OrthoHMM high sensitivity, satellite_v2 and full
OrthoFinder on four proteomes, 73,266 proteins. Do not shorten the workloads.

Periodic mode uses `measure_native_lineage_step.measure` with one-second
observations. Boundary mode uses `measure_lineage_boundary_step.measure`.
Both use the same point reader, native worker/timer, one-second completion
polling and native-step memory accounting. Both collect mandatory native
pressure, aggregate lineage counters and host CPU endpoints. Boundary mode
has exactly two observations and explicitly lacks interval screening.
This is incremental periodic-collection cost, not monitor-free or cold-cache
timing. No measured overhead is subtracted from reported tool runtimes.

Use exclusive spark-7ff0, 20 CPUs/96 GiB, one task at a time, no requeue,
fresh outputs/caches and pinned recipe/runtime hashes. Preserve the 900-second
native timeout and one-hour scheduler limit. Verify queue and host state
before release. Record every terminal controller result. Do not SSH or inspect
native outputs during the panel, stop unrelated jobs/services, or change
the frozen workloads. Freeze and validate the launcher/recipe before submission;
this plan alone does not authorize execution.

## Endpoints

For each assigned method/repetition pair report the signed quantity
`100 * (periodic_wall / boundary_wall - 1)`. Preserve the original limits:
each pair at most 10%, each method's three-pair median at most 5%, each
native command at least 60 seconds. Require all nine assigned pairs for
an overall numerical budget conclusion. Retain all 18 outcomes, including
failures, mismatches and missing measurements. Do not selectively replace
or rerun unfavorable pairs. Negative overhead is variation, not proven speedup.

Validate exact scheduler identity/resources, native commands, input order and
checksums, runtime before/after, raw/report equality and native outputs.
Compare canonical outputs within each pair. Record native CPU, native-step
memory peak and GNU-time process RSS as distinct accounting quantities.
Keep execution/provenance validity, output equivalence, overhead budgets and
environmental uncertainty as separate conclusions.

Report every periodic original/narrow flag, both modes' whole-command
screens, lineage differences and pressure diagnostics. Boundary interval flags
must remain null, not an empty list suggesting a quiet run. No flag is
dismissed based on a low overhead estimate, longer averaging window or
small signed root-minus-job difference. Lineage reads are non-atomic and
cannot identify all transient competitors or bound their effects.

## Scientific Timing Gate

This experiment does not authorize the 27 scientific scaling runs or resolve
the remaining native CPU discrepancies and during-read lifecycle behavior.
A separate prospective inclusion policy must justify environmental evidence
and measurement uncertainty before new scientific timing collection. Numerical
budget passage alone is insufficient. Retain failure and missing-data reporting
under that policy; do not retrofit eligibility to favor observed runtimes.
These three repetitions on the smallest input do not establish overhead at
larger sizes, an accuracy advantage or publication readiness.
