# Prospective Full-Node CPU Controls

Freeze before collecting DGX outcomes. The preceding native diagnostic has
21 satellite_v2 narrow flags within phylogeny and one OrthoFinder flag.
High process creation co-occurs with satellite flags but does not explain
the OrthoFinder flag. Earlier successful single-core controls did not test
this workload. No historical flag or eligibility is changed by this protocol.

## Conditions And Order

Three conditions, three repetitions each, sequentially on spark-7ff0:

1. Steady: 20 native workers, one pinned to each allocated CPU, performing
   bounded arithmetic for20seconds after a shared start signal.
2. Churn: the same20 native workers/affinities/duration, each repeatedly
   creating one child, immediately exiting that child and waiting for it
   before creating another. Count successful creation/wait pairs, require
   zero exit status for each, and retain PID and wait status on failure.
3. Contended: the same steady native workload plus one arithmetic worker in
   the observer/batch cgroup, pinned to the lowest allocated CPU, sharing
   the same20-second start signal. It is outside the target native step but
   within the authorized job, not an unrelated user's workload.

Fixed block orders: steady/churn/contended; churn/contended/steady;
contended/steady/churn. No selective retries, changes in worker count,
duration, condition order, CPU thresholds or input panel after seeing data.

Use exclusive20CPU/96GiB allocation, a60-second native-command timeout and
45-minute scheduler limit. Check the DGX queue before submission. Freeze
recipe/source hashes and pinned Python interpreter before execution. No SSH
inspection during the panel; retain terminal controller records. Do not stop
unrelated services or modify host settings. These are engineering controls,
not OrthoHMM/OrthoFinder benchmark timings or a replacement for scaling.

## Safety And Witnesses

No recursion: at most20 direct workers and20 short-lived children, plus the
native parent and one separately controlled competitor. Each churn worker
has a200,000-creation cap; reaching it invalidates the intended workload
rather than extending it. All work has monotonic deadlines. Clean up only
children created by this experiment, and wait for their termination. Preserve
failures and partial files. Local smoke tests use fewer CPUs and shorter
durations and are not counted as DGX outcomes.

Retain native-parent and every worker's PID, parent PID, cgroup membership,
allowed/pinned affinity, readiness, start/end timestamps, CPU consumption,
completed creations and exit status. Record self CPU and waited-for-child
CPU separately. The optional competitor must remain in the observer step,
on the recorded CPU, and consume5-21CPU-seconds. Do not infer achieved load
from requested worker count; report actual CPU rates and creation rates.

The common interval across all20 native workers must be at least19.5seconds.
In contended trials, competitor overlap with that common interval must be
at least19.5seconds. A dose/overlap/scope violation is a failed control, not
a negative screen response. Native work in every trial must finish within
the unchanged native collector's command boundaries.

## Observation And Analysis

Use the existing full-command dual-bracket collector and one-second cadence,
including native pressure, cgroup-frontier checks and native memory. Preserve
original and narrow screens unchanged. Retain all points from startup through
final observation, plus independent workload-ready/go/completion witnesses.
No subtraction, residual clipping or threshold adaptation.

After every trial is terminal, independently replay observations and validate
all witnesses, recipe and scheduler records. Analyze all intervals and also
the common-work subset: intervals whose enclosing original host windows lie
wholly inside the common native interval (and competitor interval for the
contended condition). Boundary intervals remain in the complete report.

For every repetition report flag counts/reasons, signed residual distribution,
native CPU rate, observed host process creations/context switches, known
worker creations, native pressure, named outside-target activity, root
residual and observer CPU. Report within-block steady/churn differences as
exploratory descriptions, not independent-interval significance tests.

For each valid contended repetition, report whether at least one common-work
interval detects excess unassigned CPU. A missing detection is retained as
failed positive-control response. For steady/churn, report every flag without
presuming an all-clean outcome or treating kernel effects as identified foreign
work. Native-only flags can test the insufficiency of a simple foreign-load
interpretation; they do not prove the cause of historical flags or validate
an arbitrary relaxed threshold.

Separate workload validity, positive-control response and screen description.
No automatic scientific timing admission follows from any outcome. Collector
overhead, larger inputs, non-CPU isolation, a frozen scientific inclusion
policy and the27-run matched scaling experiment remain separate requirements.
