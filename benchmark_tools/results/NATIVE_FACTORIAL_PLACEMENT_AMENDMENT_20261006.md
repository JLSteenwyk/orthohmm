# Native Placement Failure And Prospective Correction

## New Collector Integration Prepared

Separate native_factorial_allocated_placement.py records the initial actual
64-slot step affinity/NUMA topology, deterministic32-physical-core selection,
own-worker binding and final placement, all cgroup ancestor limits and sources.
It rejects mismatched contexts/ancestry/job/topology/cpusets/RAM or admission.
Only PID0 is bound; no scheduler policy or unrelated affinity is modified.

measure_allocated_threadripper_scaling.py launches a core-bound bootstrap,
validates ready placement before any native release, and observes the selected
OS IDs instead of literal0-31. Old collector/helpers remain unchanged. Exact
32 cores/64 slots/128GiB/85800s/1s/30s settings, abort/freshness checks, CPU and
memory scopes, host monitoring, completion and report finalization persist.
Distinct command/collector schemas and source bindings expose the amendment.
replay_allocated_threadripper_scaling.py independently reproduces selection,
ready/launch/source/command bindings and periodic affinity outcomes, and reuses
the unchanged accounting arithmetic. It requires all native/job/host/completion/
reporting evidence and refuses old schemas; it does not authorize production.

261 tests pass3.59s/zero failures/errors/skips in finalv3 XML. Earlier v1
failures expose mutable test-fixture aliasing, corrected with per-test copies;
retain v1 and passing129-test v2 stages. Scientific3.10 preparation checks all
920 original helpers unchanged. No runtime installs or frozen scientific edits.

run_allocated_threadripper_fixture.py prepares one independent engineering
route: isolated scientific3.10/four workers/six wall seconds each/4MiB each,
32-CPU inherited affinity, new collector/replay and unchanged resource endpoints.
Preparedc5e5f5b7 binds13 current source/runtime records, exact command, original
plan and safe memory/disk plus pressure/node/queue observation. Swap is full;
MemAvailable remains above128GiB; the short fixture is not a large-memory job.
It records contention rather than claiming isolation. Job23897 is submitted
once and held64CPU/128GiB/20min/no-requeue/singlebizon, not node-exclusive.
Held actual scheduler SubmitLine binds both prepared and driver57e3be59 SHA.
Initial held observer incorrectly expectsNumNodes1 instead of actual1-1;
same-job readback corrects this representation, retaining the observation error.
No duplicate/release/native inference is claimed. Push before fresh release.

Real scheduler binding and complete replay still need fixture execution and
independent terminal evidence. The full native controller/reviewer/request/
conversion/admission integration is not yet implemented. No remaining native
identity, old failure or historical measurement is restarted or reinterpreted.

## Verified Failure Disposition

The preceding goal turn makes real progress: conversion23892 completes and
assessment23894 starts. This continuation reads the active prompt/newest ledger
and rechecks same23894 RUNNING, not a stopped or blocked goal. No native retry,
new quiet-window requirement, scheduler/service change or unrelated job action.

New standalone reviewer13b37da700e3613bfc4745981d71e6d99235a9f55cab38342fd56bacc892cdf3
is outside the unchanged920-helper plan. It accepts only the exact Slurm CPU
binding refusal, actualFAILED1:0 allocation/batch andCANCELLED0:64 zero-second
launch step, original ready-timeout wrapper and abort/cleanup gates, exact
frozen command, absent native/measurement artifacts and empty host-monitor
stream. It refuses successful/running/borrowed/other failures, unexpected steps,
evidence of native execution, source/input/runtime mismatch and stale evidence.
All original limits, source bindings and failed outcomes remain unchanged.

Actual scientific Python3.10 CLI returns a distinct failure review for23891:
benchmarks/work/native_factorial_launch_failure_review_23891/review.json,
263870bytes/SHA2561f6167cf997f37dcff5f795c9479475bd065aa7359a68a73994f7f53abf4ed9a.
Runtime pre/post brackets and lookups reproduce; a fresh full runtime inventory
matches in28.247s. This is not continuous integrity or a native resource replay.
The137-second allocation duration is not inference time. Resources remainnull,
native outcome isnot_started, and accuracy/output/timing/success/automatic-retry/
publication-ready flags remainfalse. Only terminal failure resolution and the
next different identity after unchanged fresh gates are authorized.

[Independent binding](native_factorial_launch_failure_binding_23891.json)
is23558bytes/8b9292f8e8cc3e395bbf77bcd43d2f7c36156758fe492352ff8751b5b6dff2ee:
1044 unique records rechecked; source/request/plan/runtime/actualfailed state
agree. Originalreviewed_history accepts the actual ten-entry prefix including
this explicitly reviewed infrastructure failure. No index10 request/job or
incompatible-placement release is created. Earlier diagnosisdb8ce843 and its
initial cleanup-absence observation error are retained, not rewritten.

## Tested Selector

New selection helpera0f5c3a3b8804dedf2bfa11f3d0a624419e03bf3d63eda3456d84e1377409b3f
is also outside unchanged920 helpers. Given a measured64-logical-CPU step and
complete host topology, it requires exactly32 physical cores with two allocated
SMT slots each and picks the lowest OS ID from each core. It verifies membership,
uniqueness, topology/NUMA consistency and deterministic ordering, rejecting
partial, duplicate, offline, malformed or incompatible allocations. It changes
no affinity, submits nothing and authorizes no science/timing. Logical IDs0-31
are not presumed to represent32 distinct physical cores.

[Real-mask readback](native_factorial_cpu_selection_readback_20261006.json)
is16717bytes/23f490eaca28d9298ff8edf210e1c03d304cd535ba076290c590812731eb5650.
Actual failed-step mask plus twice-read current sysfs topology gives native
IDs40,49-79, mask0xfffffffe010000000000,32 distinct physical cores, NUMA node0.
Selection is unchanged when topology records are reversed. This is a retained
failed allocation/current-topology readback, not a new live allocation, binding
test or production inference. No old CPU mask or historical artifact is edited.

Final joined tests:206pass/2.28s, XML2.223s, zero failures/errors/skips.
30 selector cases and57 failure-review cases join119 existing native history/
request/reviewer cases. Synthetic complete-review tests assert every
non-admission flag and reject stale plan, runtime failure, mismatched session,
native evidence, changed request and reused destination. Test environment3.12;
actual failure review and real-mask readback use scientific3.10. No installs.
Retain earlier169/176/205-test XML as preparation-stage evidence; the final
206-test v2 XML covers the current source and added interleaved-SMT case.

## Remaining Integration

Slurm binds tasks within allocated CPU subsets; a literal OS mask cannot be
assumed available before allocation. The current official
[CPU-management guide](https://slurm.schedmd.com/cpu_management.html) and
[srun reference](https://slurm.schedmd.com/srun.html) support this design basis,
not proof of future behavior on installed Slurm24.05.2.

Before further native execution, implement and validate a prospective versioned
route that observes the actual step allocation, chooses32 allocated physical
cores and records/binds that choice. A core-bound64-slot step is a candidate
bootstrap, not an already validated launcher. Worker/native guards, thread
affinity observation, replay, reviewer, request/plan amendment and subsequent
conversion/admission bindings must agree on actual placement. Preserve original
collector/native/replay/helper files and all historical receipts; they still
require0-31 and cannot silently consume an alternative mask or be relabeled as
a successful new route. Freeze explicit source/protocol bindings before launch.

Retain unchanged scientific method/factors, input bytes/enumeration, run order,
32 physical-core native count/64 scheduler slots/128GiB, memory/CPU accounting,
capacity checks and failure handling. Record real CPU/NUMA topology and timing
comparability limits. Do not borrow unallocated cores, force node exclusivity,
disrupt other jobs, assume isolated performance or retry failed identities
automatically. This is an allocation compatibility correction, not a quiet-host
gate. The selector alone does not solve the complete production-placement issue.

Same assessment23894 remainsRUNNING at12:21 with the original six endpoints
already executing. Recheck actual terminal state and use unchanged independent
admission before reporting new scores. The full publication goal remains active
and unproven; this bounded engineering milestone is not publication readiness.
