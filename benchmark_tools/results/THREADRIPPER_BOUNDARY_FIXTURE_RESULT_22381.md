# Native Boundary Integration Result

Job 22381 was submitted once after prospective source milestone `8da47d44`
was committed and pushed. Allocation and batch completed 0:0 in nine seconds;
native step completed 0:0 in eight seconds. No restart or requeue occurred.
The command was the pinned `/usr/bin/true`, not scientific inference.
See the [submission](threadripper_boundary_fixture_submission_20260930.json),
[prospective protocol](THREADRIPPER_BOUNDARY_FIXTURE_PROTOCOL_20260930.md),
[fresh raw audit](threadripper_boundary_fixture_audit_22381.json) and
[terminal readback](threadripper_boundary_fixture_terminal_22381.json).

| Prespecified check or diagnostic measurement | Actual result |
| --- | --- |
| Native exit / timeout | Zero / false |
| Native raw points | Exactly two |
| Affinity and effective memory cap | CPUs 0-31 / 128 GiB |
| Scheduler task slots / exclusive-node slots | 64 / 192 |
| Allocation time limit | Five minutes |
| Native monotonic command interval | 0.002617555 seconds |
| Native task-subtree CPU bracket, including wrapper | 0.004793 seconds |
| Native-step lifetime peak, including launcher | 12,595,200 bytes |
| Fresh raw audit versus in-job audit | Identical |
| Rechecked direct evidence pins | 21, all matching |
| Current helper inventory pins | 832, all matching |
| Scientific timing / workload isolation admission | False / false |

The terminal controller checker verifies actual local allocation identity,
resources, command, cwd, time limit, terminal state and zero exit. Separate
accounting records confirm allocation/batch/native step status. Fresh raw
replay checks both points and command/placement/status, memory, completion,
host and finalization evidence rather than trusting the in-job result.
The [committed-source readback](threadripper_boundary_source_readback_22381.json)
separately hashes Git blob payloads from the prospective source revision: all
832 helpers and three recipe files match their exact pinned bytes/checksums.
This is source identity, not complete loaded-runtime/dependency closure.

The common host observer has two snapshots and reports competing work:
90.5822 observed foreign CPU-core equivalents over its process interval.
This is not a whole-run resource estimate or a quiet-window certificate.
The two raw snapshots retain one and zero process-read errors; zero observer
exceptions in the summary does not imply zero raw process races. No competing
job/service was stopped and no quiet-window question was repeated.

## Scope And Next Work

This completes the stated tiny native integration check of the newly added
boundary arm under this prospective recipe. It does not test descendant-heavy
inference, periodic cadence, pair-output equality, observer slowdown or full
production environmental handoff. Prior successful periodic diagnostic 22380
was not repeated; its original failures remain unchanged. No scientific
default, method, input, score or historical receipt was altered.

Native resource values above are diagnostic, not benchmark speed or algorithm
RSS. The CPU bracket contains wrapper work, the step peak contains the launcher,
and overlapping job/reporting peaks must not be added or subtracted. No total
monitor cost, runtime advantage or overhead correction follows. The five-minute
allocation does not validate the full production timeout/reporting reserve.

The 54-task paired engineering plan and 27 production timing identities remain
unstarted. They still need a complete independent pair audit, prospective stable
source/runtime and environmental execution recipe, applicable overhead checks
and a verified quiet local window. Do not silently upgrade the historical
plan's pending implementation binding: a future recipe must explicitly bind
this implementation and the retained integration result. Other QfO uncertainty,
source/rights and release requirements remain open. Full publication goal active.
