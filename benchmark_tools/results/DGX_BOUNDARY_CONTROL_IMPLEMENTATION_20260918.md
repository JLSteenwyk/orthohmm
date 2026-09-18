# Boundary-Only Control Implementation

`measure_frontier_boundary_step.py` is a separate control, not a change to
any frozen collector. It retains the frontier collector's same native
worker, scheduler/resource checks, start handshake, one-second completion
polls, timeout/cleanup, pre-command and post-command frontier observations,
and final step-memory observation. It omits periodic counter reads and their
point-file serialization while the command is running.

This enables a prospective estimate of incremental periodic observation
cost. It does not represent an uninstrumented command: worker snapshots,
polling, start/end counter reads, and wrappers remain in both arms. In
particular, completion-detection latency is outside the worker's retained
native command-wall interval in both arms.

Boundary reports validate complete-command enclosure, unchanged whole-span
CPU thresholds, stable scope identities and monotonic counters. They set
`interval_screening_available=false` and `flagged_intervals=null`, never
an empty list implying passed interval checks. Periodic screening cannot be
inferred from the two endpoints. Scientific timing admission remains false.

Tests exercise mechanical worker/lifecycle equivalence to the frozen
collector, omission of reads across multiple completion polls, native
success/nonzero/timeout propagation, invalid boundary evidence, and cleanup
after a read failure. Replaying first/last points from the three retained
native runs checks boundary arithmetic, but does not turn those instrumented
runs into control-arm observations. No paired native overhead result exists
yet for this collector.

Before execution, freeze a complete paired panel, command/input/runtime
identities, balanced ordering, repetitions, engineering budgets and failure
handling. Include representative native input sizes; the existing645-protein
fixture alone cannot establish overhead for the scaling workloads. Native
timeouts in this control remain capped at900seconds by the original worker
validation. Longer commands require a separately reviewed prospective change,
not a silent extension or reuse of incomplete measurements. No outcome-driven
repeats, retrospective policy relaxation, or overhead subtraction is allowed.
