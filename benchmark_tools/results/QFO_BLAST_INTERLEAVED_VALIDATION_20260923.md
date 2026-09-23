# Interleaved Recovery Validation Amendment

## Scheduling Issue

The [frozen search-recovery protocol](QFO_BLAST_RECOVERY_PROTOCOL_20260923.md)
remains unchanged for query membership, full database, native command,
runtime, resource allocation and output validation. The original admission
array additionally waited for the entire search array. Releasing only task
0's whole-array dependency made it eligible but did not ensure a resource
window: a 900-GiB search plus a 96-GiB scoring job leaves insufficient
scheduler-accounted memory for a 64-GiB validator on the 1,030,000-MiB node.
This observation concerns reservations, not actual RSS or causal CPU effects.

## Existing-Job Dependency Amendment

Applied only to pending tasks in the existing recovery arrays. Native
22103_0 is already complete; 22103_1 was running and was not modified.
No job was canceled, resubmitted, requeued or restarted. The native and
admission executors and original submission scripts remain unchanged.

- Admission 22105_0 remains eligible after its independently confirmed
  completed native task, retaining its in-process admission checks.
- Admission 22105_1 requires both successful native 22103_1 and successful
  admission 22105_0.
- Admission 22105_i for i=2..19 requires successful native 22103_i.
- Native 22103_i for i=2..19 requires successful admission 22105_(i-1).

Thus no third or later native batch can start until preceding completed
work has passed its independent validator. The original serial array limits
remain one; there is no increased native concurrency. A failed validator
blocks its downstream searches for review instead of being automatically
retried or bypassed. Unrelated jobs and machine configuration are untouched.

The [scheduler receipt](blast_recovery_interleaved_dependencies_20260923.json)
retains 36 bulk updates with before/after records, the separately recorded
task-1 admission amendments, and final snapshots of all 38 pending tasks.
After every bulk update, array identity, throttle, CPUs, memory, time limit,
requeue/restart counters, command and working directory were unchanged.
An independent traversal checked the final exact dependency sets and found
no cycles. All scheduler updates returned success.

## Scientific Boundaries

Dependency changes create opportunities for validation; they do not prove
that any validator ran or passed. Final batch/whole-panel admission, retained
prefix review, ordered merge, conversion, clustering and scoring remain
required. Original interrupted output is immutable and the old downstream
chain remains held. These recovery costs are not matched comparative timing.
No DGX activity or change to the scientific parameter neighborhood is involved.
