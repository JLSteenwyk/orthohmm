# Prospective Bracketed CPU Screen

## Purpose And Scope

This freezes an engineering diagnostic before any new timing panel. It does
not authorize another27-run panel or establish scientific timing admission.
The earlier native smokes verified functional operation but their counter
windows cannot support host-minus-native arithmetic.

Required order for a prospective measurement:

1. Start a distinct native Slurm step and publish its identity, but do not
   take the native pre-work counter snapshot yet.
2. Finish the observer's host pre-snapshot, then release the native worker.
3. Take the native pre-snapshot, execute the unchanged command, and take the
   native post-snapshot before announcing completion.
4. Take the observer's host post-snapshot after native completion. Preserve
   the step until all required final counter evidence has been acquired.

Require complete monotonic read brackets, stable boot/topology, separate batch
and native scopes in the same job, no counter-read errors, nondecreasing CPU
counters, and the command fully enclosed by native reads. The host window
must fully enclose native reads. Reject unavailable or mismatched evidence;
do not infer missing intervals from neighboring observations.

## Frozen Diagnostic Arithmetic

Compute signed unassigned CPU as outer-window host busy CPU minus inner-window
native cgroup CPU. Host busy sums user,nice,system,irq,softirq in SC_CLK_TCK
units without double-counting guest. Cgroup CPU uses integer usage_usec
differencing before conversion to seconds. Do not subtract observer CPU:
including it makes the residual operationally conservative, but does not
turn it into a mathematically proven bound.

Flag unassigned CPU above0.25average cores, using native command wall time
as denominator. Flag signed residual below-0.5CPU seconds as accounting
discrepancy. Flag any host steal increment. Preserve the signed residual;
do not clamp negatives, rescale native time, or claim statistical significance.
The0.25-core and0.5-second constants are engineering screens, not calibrated
confidence limits. Passing leaves controlled_workload_verified and
scientific_timings_admitted false.

## Actual Historical Check

Executed `screen_bracketed_cpu.audit_existing` on the retained three native
smokes. All were rejected for missing outer bracketing, as expected. Observer
pre-read finish minus native pre-read start was7.534166ms for high-sensitivity,
7.471046ms for satellite_v2, and16.295519ms for OrthoFinder. These are measured
window overlaps, not estimates of foreign CPU. Exact snapshots and source-file
checksums are in `dgx_native_window_rejections_20260918.json`.

The three completed native jobs remain valid functional smokes. Neither they
nor the original27timing runs receive new scientific admission through this
audit. No scientific run was restarted or changed.

## Remaining Gates

Implement the prospective handshake and test quiet, known short-lived and
sustained sibling CPU controls with raw evidence before using the diagnostic.
A run-wide average can conceal concentrated bursts, and CPU counters cannot
exclude memory, I/O, thermal or frequency interference. Interval-level
observation and an explicit treatment of accounting/read-window uncertainty
remain necessary for a new scientific inclusion protocol. That protocol must
also freeze resources, memory scope, command boundaries, run order, repeat
policy, missing-run treatment and reporting before replacement times are
examined. No claim of rigorous interference bounds is made here.
