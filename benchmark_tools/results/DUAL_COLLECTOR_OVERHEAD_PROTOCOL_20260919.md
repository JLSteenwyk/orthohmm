# Prospective Dual-Collector Overhead Panel

Freeze before collecting new paired timings. The full-node controls validate
the intended known-competitor detection response, not monitor overhead or
scientific timing eligibility. Historical native residual flags remain.

## Design

Run a complete fresh18-task panel: three methods, three repetitions, two
arms per pair. Derive method commands, four-proteome inputs, seeds, native
output semantics, resource limits and exact task order from the frozen
pressure-overhead-v2 plan SHA-256
`b644e165dbf4d0beabf1cf4d9b6c314de522e3ebd1b91598ebebea99094c8fff`.
Relocate only output/cache paths. Methods are frozen OrthoHMM high sensitivity,
OrthoHMM satellite_v2 and OrthoFinder full. Do not substitute shorter workloads.

The periodic arm uses `measure_native_dual_bracket_step.measure` at the same
one-second cadence. The boundary arm uses the existing
`measure_frontier_boundary_step.measure(native_pressure=True)`. Both retain
the same command worker, native timer, one-second completion polling, native
memory and pressure endpoints. Boundary collection has no periodic frontier
scans. Dual validation done before command release and after completion is
outside the native timer; this comparison estimates incremental periodic
collection cost, not zero-monitor or cold-cache runtime.

Use exclusive spark-7ff0 allocations,20CPUs/96GiB, no requeue, sequential
execution, pinned Python and runtime/recipe hashes, fresh output/cache paths
and retained terminal controller records for all tasks. Check queue and host
activity before starting. No SSH during the panel; no unrelated jobs or host
services may be stopped to obtain favorable results. Existing frozen native
timeout and input order remain unchanged. Capture failures as failures.

## Endpoints And Gates

For every assigned pair retain native wall durations and report
`100 * (periodic_wall / boundary_wall - 1)` with its sign. Preserve the prior
numerical budgets: each pair at most10% overhead and each method's median
at most5%. Both runs must take at least60seconds. Require all nine assigned
pairs; do not replace failures with selected historical runs or repeat only
unfavorable pairs. A negative difference is variation, not proof of speedup.

Bind scheduler identity/resources, exact native argv, input checksums/order,
runtime before/after and raw report files. Validate outputs and same-method
canonical output equivalence within each pair. Record actual CPU and consistent
native-step memory accounting; do not conflate process RSS with cgroup peak.

Keep four conclusions separate: execution/provenance validity, output
equivalence, numerical overhead-budget result, and environmental uncertainty.
Report every original and narrow periodic flag and both arms' whole-command
screens, native pressure and frontier diagnostics. No flag is deleted because
overhead meets its budget. Boundary arms cannot exclude brief interference.
Numerical budget passage alone cannot establish clean comparative timings.

Do not use longer-window descriptions to redefine an accepted run after
seeing this panel. Any scientific inclusion policy must be frozen separately
before the27 matched scaling runs, justify how observed foreign activity and
non-atomic accounting uncertainty are treated, and retain failed/missing data.
This overhead experiment does not itself authorize those scientific runs.

## Reporting

Retain all18 outcomes,9 pair identities, signed differences, per-method
medians, native output counts/canonical identities, runtime checks, raw CPU
screens, memory/pressure, scheduler records and complete failure reasons.
Report any missing evidence as incomplete, not clean. Keep all preceding
panels and their failures unchanged. Do not subtract measured overhead from
published tool runtimes or claim accuracy/general superiority from this panel.
