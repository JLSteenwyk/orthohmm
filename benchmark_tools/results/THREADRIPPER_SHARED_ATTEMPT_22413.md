# OrthoFinder Repeat Aborted Before Native Inference

This retains the original failure review. The subsequent
[repair and explicit history resolution](THREADRIPPER_PRENATIVE_REPAIR_20261004.md)
does not change any failed verdict or supply missing native endpoints.

Index 17, full OrthoFinder 3.1.5 on twelve proteomes, repeat 1, was released as
22413 after the canonical index-16 review passed. Fresh accounting reports
parent and batch FAILED 1:0 at 3:03; the parked native step completed 0:0 at
0:36. The parked step's zero exit is not successful OrthoFinder inference:
its gate received `{"abort": true}` and wrote `observer_did_not_release_native`.
There is no native log, completion record, pressure-point series or OrthoFinder
output directory. Native wall/CPU/peak-memory endpoints are absent, not zero.

The [machine-readable failure review](threadripper_shared_prenative_failure_22413.json)
is 9,252 bytes, SHA-256
`9d643efb649b5a7c67cdda51fc6e359b5baaa89d938fce5e45c276d3186ea115`.
The [review source](review_shared_prenative_failure_20261004.py) binds the exact
request, plan, preflight, worker lifecycle, wrapper failure and retained process
stream. It checks abort/no-native-artifact consistency and retained pre/post
runtime-lookup decisions without rerunning full runtime checks or inference.
The [successful terminal controller query](threadripper_terminal_controller_22413.json)
was observed before purge and transcribed exactly from the tool output;
observation timestamps were not captured and are not invented. The review
corroborates its state/exit/start/end with a fresh `sacct` query.

## Failure Mechanism

The worker's retained error is `Frozen input/source identity changed`, naming
`measurement/host_processes.jsonl`. The collector begins its periodic observer
before release checks so those checks cannot delay its first sample. The
preflight helper, however, hashes the entire append-only stream as immutable
evidence and later rechecks that hash. Its second collector check also requires
the stream to still have exactly one sample. These assumptions conflict with
legitimate observation progress during a slower release path.

Retained collector sample starts are 738965.782785931 and 738995.783825298
monotonic seconds, separated by 30.0010393669 seconds. The worker's two
preflight snapshots span 738992.394661061 through 739001.555586249: the periodic
append overlaps preflight. The hash check fails before native release.
This supports a mutable-stream race, not an OrthoFinder algorithm failure,
rejection of ordinary background CPU demand, or insufficient memory capacity.
The worker records available memory of 332,716,515,328 and 332,481,789,952 bytes,
both above the frozen 137,438,953,472-byte minimum. Observed foreign demand is
42.0544714236 CPU-core equivalents; contention remains accepted diagnostically.

The wrapper records the release-check ValueError. Context cleanup subsequently
raises `Environmental worker was not successfully joined at release`, masking
the earlier error in the enclosing executor log. Both errors and the child's
terminal exit 1 remain retained. No passing preflight is fabricated.

## Next Actions

Seventeen focused tests verify this classification, refuse borrowed/contradictory
identities, gates, artifacts, runtime decisions and capacity, and reproduce the
current helper's hash failure after a normal second stream append. They document
the defect; they do not implement or validate its repair.

Preserve this failed attempt and its request/output paths without retry or
overwrite. Repair the preflight to bind an immutable initial observation while
validating the live observer/parked worker and accepting legitimate later
samples. Preserve cadence, fresh-capacity, runtime and release safeguards;
do not move the observer after release or loosen monitoring bounds. Test the
repair, refresh affected source/runtime/protocol bindings and explicitly resolve
history before index 18 (full OrthoFinder, four proteomes, repeat 2). The current
continuation/reporting helpers do not support this new abort outcome: do not
forge their expected successful-native or cadence-only resolution schemas.

The [current resource snapshot](threadripper_shared_panel_snapshot_20261003_v17/panel.json)
contains seventeen resource-reviewed attempts, not this eighteenth pre-native
abort. Its missing index-17 endpoints must remain distinguished from unrun work
in final reporting. No quiet-window requirement, unrelated workload change,
scientific retuning, native accuracy claim or publication-readiness claim occurs.
