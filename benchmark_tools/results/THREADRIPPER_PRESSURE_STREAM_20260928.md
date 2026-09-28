# Pressure Review Throughout Native Execution

The bound post-measurement review now includes the existing per-point CPU,
memory and I/O PSI observations. No new native sampling, scientific settings,
input bytes or run order were introduced. The executor requires explicit
`maximum_pressure_percent` limits for all three resources and a positive
`maximum_pressure_sample_period_s` before creating the attempt. These limits
are not automatically inferred or approved from historical measurements.

The pressure evaluator verifies both host brackets at each point, including
read errors, boot, CPU set, observer membership, monotonic observation windows,
PSI syntax and nondecreasing counters. It uses the second host bracket of each
point for successive cadence-scale rate comparisons, not the much shorter
within-point intervals. Every such interval is checked: an intermediate burst
cannot disappear into an end-to-end average. Point duration and sampling-period
bounds are explicit, and the observations must bracket native execution.

The file-backed audit requires contiguous numbered point files, binds each
file's bytes and checksums, rechecks them after evaluation, and verifies that
the file inventory did not change. Only individual decoded points are streamed;
the path/checksum manifest still grows with point count, as does the retained
audit report. Additional post-native reads and hashing require full-scale
budget/overhead validation. Previous fixtures do not establish that budget.

The final report separates `sampled_process_policy_satisfied` from the nested
pressure verdict and their conjunction `sampled_environment_policy_satisfied`.
The executor checks the conjunction, retains negative reports and fails the
attempt without retry or timing admission. The original standalone process
status retains its narrower meaning.

## Verification

262 focused pressure, process-stream, collector, executor, replay, wrapper and
lifecycle tests pass. Pressure controls include bursts in each resource,
missing/corrupt reads, counter decreases, changed boot/CPU set/scope, overlapping
observations, gaps and absent native bracketing. A file-backed control verifies
that clean process observations cannot override failed pressure evidence.

The [historical negative replay](threadripper_pressure_stream_negative_replay_20260928.json)
checks all raw point hashes for the retained high-sensitivity, phylogenetic and
OrthoFinder fixtures. Deliberately zero-percent PSI limits reject all three.
These are negative-control settings, not prospective production policy. Maximum
cadence-scale CPU `some` midpoint percentages were 20.7069%, 6.2230% and 0.3452%;
I/O percentages were 0.01811%, 0.02006% and 0.02053%. These are observations from
different historical shared-host windows, not comparable method performance.

System pressure includes the benchmark's own work. A failed bound does not
identify foreign interference, and zero recorded stalls do not prove isolation.
Midpoint rates are approximate because reads are non-atomic and accounting can
be delayed; values are not clipped. System CPU `full` is not interpreted.
The historical PSI parser was reused without altering DGX collectors or accessing
the DGX. See its retained [Linux PSI documentation](https://docs.kernel.org/accounting/psi.html).

Whole-run executable/configuration stability, reviewed prospective limits,
updated source freeze, native/full-scale validation and a quiet window remain
open. No production benchmark or unrelated service change occurred.
