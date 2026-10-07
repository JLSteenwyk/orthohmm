# Allocated Native Review Handoff

## Prospective Postprocessing

Native identity10/P1C0R1 is already running as23902. No allocated terminal
review is queued at this inspection. Add only a small scheduled batch wrapper
around the existing frozen `review_allocated_native_factorial_attempt` CLI;
do not modify any of the 14 bound allocation-route sources, historical plan,
920 helpers or native scientific method. The original reviewer retains all
terminal/runtime/resource/environment/scientific-output gates.

[New batch](allocated_native_review_23902_20261006.sh) requests two CPU slots,
128GiB, six hours, one task/node onbizon and no requeue. The larger memory
is review-only, consistent with the successful retained old full-review
postprocessing envelope; it does not change native23902's limits. Reject
unscheduled, wrong-CPU or same-native-job invocation. Sanitize Python/loader
variables, preserve the original scientific Python3.10 venv entry point,
enable unbuffered/faulthandler diagnostics, and log/check available RAM >=128GiB
at execution. Unsafe capacity stops this attempt; no automatic retry or
unrelated-work intervention is authorized.

Fixed request1355ae3b binds original native23902 and identity10. The fresh
review destination is `benchmarks/work/allocated_native_factorial_terminal_review_23902_v1`.
The frozen reviewer itself requires actual terminal accounting before creating
its destination. It reviews retained failed outcomes as such where possible;
process success alone is not accuracy or resource admission. A failed/partial
review, signals, timeout, diagnostics and scheduler outcome must be retained.
No historical failed reviewer or native identity is overwritten or repeated.

Commit/push this batch before submission. Submit one held job with
`afterany:23902`, fresh log/destination and no array. Inspect its actual
owner, held state, dependency, command/cwd, two CPUs/128GiB/six-hour envelope,
no-requeue/zero-restarts and source/input hashes before release. Pin the held
submission in an immutable receipt and scheduler comment. Release that same
handle once after the fresh checks, then observe it as dependency-pending,
not completed review. Do not submit another reviewer on continuation.

## Validation

[24 focused tests](allocated_native_review_batch_tests_20261006_v1.xml)
pass in1.02s, zero failures/errors/skips: seven new batch checks and the
existing17 allocated-review tests. Bash syntax and real CLI-help loading in
the scientific Python3.10 environment pass. These are scheduling/component
checks, not an actual full postprocessing review or missing native score.
Before any submission, strengthen the RAM gate to an explicit `require`, not
an optimization-dependent Python assertion. Retained v1 tests precede this
change; [v2 validation](allocated_native_review_batch_tests_20261006_v2.xml)
passes the same 24 cases in0.90s with no failures/errors/skips and a static
guard against reintroducing `assert n`. No scientific reviewer is changed.

The current native job remains live and will not be restarted. Conversion,
assessment and independent accuracy admission require its actual reviewed
outcome, not this scheduling milestone.11/12 remain sequential behind valid
reviewed history; retain9failure without automatic retry. The full publication
goal remains active and incomplete.

All eventual timing observations are shared-Threadripper measurements under
competing CPU, memory-bandwidth and I/O demand, with unknown, potentially
tool-dependent effects. Review is separate from inference/conversion/scoring;
no isolated-performance ranking, quiet/DGX requirement or readiness claim.
