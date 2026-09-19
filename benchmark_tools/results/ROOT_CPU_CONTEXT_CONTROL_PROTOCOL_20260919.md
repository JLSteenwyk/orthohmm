# Prospective Root CPU Context Controls

## Question

Do the new root/system.slice/user.slice/init.scope observations distinguish
a known owned user-service CPU workload from the same native-only controls?
The completed lineage panel localized signed counter differences near the
root-to-system.slice boundary, but did not establish a process-level cause.
The completed read-crossing control demonstrated why partial-span signed
differences are not atomic interference bounds. No historical flag or
scientific timing eligibility changes under this protocol.

## Fixed Panel

Use four conditions, three repetitions each, in this fixed order:

1. Idle, steady, churn, user-contended.
2. Steady, user-contended, idle, churn.
3. Churn, idle, user-contended, steady.

Idle has a sleeping native child for 20 seconds. Steady and churn reuse
the bounded full-node workload: 20 workers pinned one per allocated CPU,
20 seconds after a shared go signal, and the existing 200,000-creation cap
per churn worker. Preserve readiness, actual CPU, creations, exit status,
membership and affinity witnesses. Existing workload validity checks remain.

User-contended uses the same 20 steady native workers plus exactly one
owned arithmetic worker in a uniquely named user-manager service, pinned
to the lowest allocated CPU and sharing the start signal. Reuse the existing
bounded competitor payload; do not create an unbounded burner. Require
5-21 process-CPU seconds and at least 19.5 seconds of overlap with the common
native worker interval. Record membership outside the Slurm job under the
observed user-manager subtree, and unchanged one-CPU affinity. Apply
RuntimeMaxSec=45s, MemoryMax=256M, TasksMax=8 and --collect to the owned unit.
Retain service exit and disappearance evidence. Never stop unrelated services.

Reserve all 20 DGX CPUs exclusively with 96 GiB memory and a 15-minute
outer job limit. Each native command has a 60-second timeout. Verify committed
source/interpreter identities and an empty DGX queue before submission. No
SSH inspection or transfers during the complete panel. Collect all outcomes
after termination, including failures, without selective replacement.

## Observation And Evaluation

Retain the unchanged lineage measurement and its original/narrow screens.
Add separate root-context observations at the same one-second cadence,
preserving their actual timestamps and durations. The supplementary reads
must not be represented as belonging to the existing, narrower native read
windows. Require stable boot and named-scope identity; missing counters and
decreasing counters are explicit failures with partial evidence. Root PID
membership churn is retained, not interpreted as a counter failure.

Record root/system/user/init cpu.stat, enclosing host CPU categories, clock
tick frequency, and root membership before/after each context snapshot.
The three named children are not an exhaustive root-child inventory. PID
snapshots do not identify per-task CPU, PID reuse or unobserved transient
tasks. Signed residuals remain signed; no clamping or timing correction.

After all 12 conditions finish, independently replay raw observations and
validate all workload/overlap/service/source/scheduler witnesses. For each
valid user-contended repetition, require at least 5 CPU seconds in the
user.slice aggregate over the enclosing workload span. A failure remains a
failed positive-control response, not a relaxed threshold. This response
does not establish exact recovery of the known load or exclude other user
activity. Report all idle/native-only outcomes without presuming zero flags.

Report all-interval and wholly common-work-interval distributions separately:
named-scope CPU, signed root-minus-system and root-minus-three-child values,
enclosing host categories, PID membership changes and original/narrow flags.
Do not sum overlapping host guest/user categories or subtract observations
with different windows to claim causal decomposition. Describe fixed-block
condition differences without treating dependent intervals as independent
replicates. Preserve unknown coverage and all failed trials.

These are observer engineering controls, not native accuracy or speed
benchmarks. Context-enabled overhead, actual native-tool attribution,
non-CPU isolation and prospective scientific timing inclusion remain
separate. No scientific run is admitted automatically by any outcome.

## Current Implementation Status

`probe_root_cpu_context.py` provides the read-only observation and comparison
primitive. Its 63 focused tests cover scope identity/order, signed residuals,
host counters, membership churn and partial failures. The complete panel
runner, integration and independent workload replay are not yet implemented
or deployed; no DGX result is claimed by this protocol freeze.
