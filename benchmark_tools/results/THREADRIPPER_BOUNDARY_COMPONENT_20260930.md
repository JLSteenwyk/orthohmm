# Native Boundary Control Component

Implemented the boundary arm required by the
[prospective paired-control plan](THREADRIPPER_NATIVE_OVERHEAD_PREPARATION_20260930.md).
This is component-level engineering progress, not an executed overhead experiment.
The user's quiet window remains deferred; no scheduler submission, new host-load
probe, DGX action, native inference or interruption of unrelated work occurred.
An empty current Slurm queue does not establish a quiet host.

## Implementation

[Boundary collector](../measure_threadripper_boundary.py) exposes a callable
measurement function, not a launch CLI. It reuses the unchanged periodic
collector's native worker, command clock, timeout, placement and release gate.
Frozen settings remain 32 affinity CPUs, 64 scheduler slots, 128 GiB, 85,800-second
native timeout and one-second completion polling. It collects exactly two
native points: before command release and after completion. Both arms retain
the common 30-second asynchronous whole-host observer.

The collector retains raw points, command/placement, native completion, step
and job memory, host observations and reporting-finalization evidence. Collection
exceptions retain a failure marker when the guarded collection phase is reached.
Absent or incomplete receipts cannot pass the checker. Native nonzero exits
and timeouts remain failures, even when their resource receipts can be replayed.
All scientific timing, controlled-workload, output-validation and publication
admission flags remain false. No production executor or plan authorization changed.

[Raw checker](../replay_threadripper_boundary.py) requires explicit external
command, launcher and worker bindings. It reconstructs screening/context,
anchor-only boundary completion and host summary from raw evidence; checks
two-point inventory and hashes, resource identity, native status/clock windows,
memory scopes/events/monotonic peaks and reporting finalization; and rejects
missing, contradictory, indirect or mid-replay-changed evidence. Relocated
archives require an explicit original directory binding. It shares low-level
arithmetic helpers with the existing checker, so this is not an independently
implemented resource estimator or a complete source/runtime/scheduler audit.

## Executed Validation

195 focused tests passed in 3.68 seconds under the existing shared Python
3.10 environment. There are 66 new tests across the two new test modules.
Coverage includes exactly two native points regardless of completion-poll count,
success/nonzero/timeout status, denied or stale release, collection failures,
surviving descendants, wrapper failure, settings rejected before directory
creation, raw-record tampering, relocation and mutations during replay.
Three composition tests feed collector-generated receipts to the actual checker
for zero and nonzero outcomes, retaining actual low-level evaluation functions.
All native launch, clocks and host/process fixtures in these tests are synthetic.

During test development, failures exposed two fixture mistakes: task-level memory
scope instead of step memory scope, and a mocked ready read without a ready file.
Fixed the fixtures, not the checker's acceptance rules. No native attempt failed
or was retried, and no scientific result changed. The existing periodic collector
and raw checker remain unchanged.

```sh
python -B -m pytest -q \
  tests/unit/test_measure_threadripper_boundary.py \
  tests/unit/test_replay_threadripper_boundary.py \
  tests/unit/test_measure_threadripper_scaling.py \
  tests/unit/test_replay_threadripper_scaling.py \
  tests/unit/test_native_completion.py \
  tests/unit/test_report_finalization.py \
  tests/unit/test_threadripper_job_memory.py \
  tests/unit/test_periodic_host_observer.py \
  tests/unit/test_prepare_threadripper_overhead.py
```

## Remaining Gates

A prospective, source-bound native integration check is still required for this
new arm. Do not relabel the historical plan's pending binding or successful
periodic calibration 22380 as boundary-arm validation. Full pair auditing must
check output equality, resources, terminal provenance, all attempts and complete
cells before applying the prospectively defined engineering budget. No paired
workload or slowdown result exists here. The contrast will not isolate the common
host monitor's cost, and two native points cannot establish continuous affinity,
containment, host isolation or environmental policy approval. Memory peaks
overlap and must not be added, subtracted or used for timing correction.

Environmental handoff, stable execution recipes, the quiet window, controlled
production timing and other publication requirements remain open. No quiet-window
question needs repeating until there is a new scheduling need. The full goal
remains active.
