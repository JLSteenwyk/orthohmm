# Corrected Replay Dispatch

## Handoff

The corrected HMM inference21706_0, native admission21720 and replay-plan
preparation21722 are already queued. `launch_qfo_corrected_replay.py` adds
the dependency-aware handoff to the existing checked replay runner. It does
not change frozen scientific settings, clustering interception, worker code
or the independent replay auditor.

The handoff requires terminal successful2-CPU/64-GiB preparation and native
admission jobs onbizon before reading the plan. It verifies clean frozen
preparation/replay executor188fde21860a70da55fee1358177485c970a11a3 and HMM
admitter7b5214a5d2169338c4cbe80293fae491ee5b951f, the exact plan/admission/
primary source identities, admitted checkpoint/input scope, fixed paths,
Python identity, command parameters, environment and four expected stages.
`corrected_evidence` rechecks the input and checkpoint identities. The
existing runner still performs its complete scientific runtime verification
and four checked clustering boundaries before producing usable output.

Only the frozen `run_qfo_corrected_replay.py` is executed. The unwrapped
`native_command` is never dispatched directly. The original plan retains
`execution_authorized: false`; a separate
`benchmarks/work/qfo_corrected_replay_dispatch_20260918.json` records the
validated dispatch, not completed execution. The wrapper uses `execv`, so
the existing runner retains the scheduled job identity and exit status.
Existing replay output or a previous dispatch record is not overwritten.

Frozen settings remainBLOSUM62, CPM0.1, Leiden seed4,32CPUs, one profile
iteration and minimum one species. Output remains
`benchmarks/results/qfo_corrected_checked_replay_v1`. Native/replayed
partitions must be compared; nonequivalent partitions cannot share scores.
Cached replay timing onbizon remains incremental/shared-host, not matched
end-to-end scaling evidence.

## Verification

30 new tests cover source/path drift, scientific settings, interpreter,
job resources, existing output, revision checks, check-only behavior and
dispatch to the unchanged checked runner. Together with existing preparation,
runner, worker-evidence, stage audit, admission and batch tests:132 passed
in2.61s. The actual pending21722 preflight rejects before reading a plan.
The current interpreter equals the one in the frozen primary command plan.

The new dispatch batch requests32CPUs/192GiB/24h/bizon, no requeue, and must
depend on afterok:21722. Independent replay admission uses the already
frozen09dec90118e9280990295f8aa9c1aa1a9171714e at
`benchmarks/work/publication_qfo_corrected_replay_admission_v1`, with a new
2-CPU/64-GiB/24h batch depending afterany of replay. It verifies terminal
execution before accepting parent/worker reports, stage partitions, input
coverage or native/replay agreement. No replay score is produced by either
dispatch or admission; candidate preparation and downstream scoring remain.

Submission identities are recorded after freezing and queue confirmation.

## Submitted

- Replay **21756**, submitted scheduler time `2026-09-18T10:10:17`,
  afterok:21722,32CPUs/192GiB/24h/bizon, no requeue.
  Dispatch executor7df7aae0a9c69ce88e28fb4fc8ae00d63cfe1a41 at
  `benchmarks/work/publication_qfo_corrected_replay_dispatch_v1`.
  The runner remains188fde21860a70da55fee1358177485c970a11a3.
- Admission **21757**, submitted scheduler time `2026-09-18T10:10:28`,
  afterany:21756,2CPUs/64GiB/24h/bizon, no requeue, using the unchanged
 09dec90118e9280990295f8aa9c1aa1a9171714e auditor.

Both resource requests and dependencies were confirmed by `scontrol`.
The actual frozen auditor rejects pending21756 at its terminal scheduler
gate before accessing partial outputs. Its future output is
`benchmarks/work/qfo_corrected_replay_admission_20260918.json`.
No corrected replay, partition-equivalence result, candidate arm or score
has completed at this submission milestone.
