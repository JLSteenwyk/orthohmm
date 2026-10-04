# Preparation Synchronization Repair And Resolved History

This follows the [retained index-20 failure](THREADRIPPER_SHARED_ATTEMPT_22416.md).
Fresh accounting confirms 22416 is FAILED 1:0, not live inference; the queue
contains no owned timing job. The preceding turn is progress (commits
`7de68449` and `0ad646b1`), not an unresolved user-input or quiet-window blocker.
No native attempt is retried. Original failed receipts and earlier reporting
artifacts remain byte-identical.

## Prospective Software Repair

The observer now publishes `environment_worker_prepared.json` only after
deployment/history/policy/source checks finish. The owner waits for this
atomic receipt before starting the measurement wrapper and parking a native
worker. It checks request/policy hashes, index/job, owned PID, boot, timestamp,
live allocation membership and stable receipt bytes. Preparation is not
preflight admission and cannot authorize native release.

The actual memory/process/pressure observations still follow the release
request. The release response remains 20 seconds, the parked-worker gate
45 seconds, and the existing preparation wait 4,000 seconds. Native timeout,
CPU/RAM enforcement, budget/freshness, periodic sampling, immutable first
sample and failure cleanup are unchanged. Preparation failures start no
collector; only the owned child is reaped. No unrelated jobs/services change.

Exactly five top-level helpers change: observer, owner, executor, panel
progress and history binder. History validation accepts only two explicit
pre-native failure kinds, each requiring its matching failure audit and
successful source-bound repair. Deadline resolution additionally requires
the original expired-preflight flag and response timeout. Both preserve
FAILED 1:0, `not_started`, failed environment/unresolved resources and no
retry, eligibility or native endpoint admission. Other failures are not waived.

Initial focused runs pass 291 cases in 24.72s and 453 in 30.46s. The
[captured joint regression](threadripper_preparation_sync_validation_20261004.json)
passes **655 tests in 35.91s**, zero failures/errors/skips, with source/test/log/
JUnit hashes checked before and after execution. It is 8,829 bytes, SHA-256
`dbe9f01229a716933e35e8c99968deadff36215456ad3253dde300cf7b849e40`.
The [immutable-sample revalidation](threadripper_immutable_sample_revalidation_20261004.json)
uses that same actual execution, not a second run. Tests include slow setup,
owned-child readiness/cleanup, malformed or borrowed receipts, no collector
after preparation failure, actual deadline classification, both resolution
contracts and the retained sample/cadence/resource/controller controls. This
is software regression, not production handoff or causal overhead evidence.

## Actual History And Runtime

The [new history result](threadripper_preparation_sync_history_20261004.json)
checks the real 21-identity prefix and selects index **21**: high-sensitivity
OrthoHMM, eight proteomes, repeat 2. It is 208,239 bytes, SHA-256
`50663773217c2ddf257a3219f04f2554f34eeb3ed48765e0cdaecda86be557f8`.
New index-17 resolution uses current immutable-sample regression while
retaining its original session, audit and decisions; its old resolution is
linked and preserved. New index-20 resolution retains the separately audited
deadline abort and fresh terminal corroboration. Neither supplies resources.
All-prefix binding passes; no successful index-20 summary is invented.

The [refreshed lookup](threadripper_private_lookup_preparation_sync_20261004.json)
is 6,661 bytes, SHA-256
`7e3ac02b46cf1ead426ec19889fc5bd23e6a0e01dbd3737bf374c3ae3180ed1d`.
Exactly those five tested helper identities differ from the prior native/OS
manifest. Scientific/private runtime, baseline, controller and command-plan
identities are unchanged. All **57,959** records match before/after startup;
native import comparisons retain **913 OrthoHMM / 1,563 OrthoFinder** modules.
This refresh is justified by concrete helper changes; reuse it rather than
repeat startup/calibration on resumption. It does not establish a live handoff.

## Reporting Checkpoint

New dated reporter/plotter preserve their earlier hash-bound counterparts.
[Snapshot v21](threadripper_shared_panel_snapshot_20261004_v21/panel.json) is
290,517 bytes, SHA-256
`9ded273a1383ca4cf027ed40bc28310f553cad287f24f1ce0892784cd3425af0`:
21 reviewed attempts, 19 measured, 18 eligible, exclusions `[0, 17, 20]`,
two explicit pre-native aborts and six not-yet-reviewed identities. All first
twenty rows are unchanged. Both four-proteome OrthoHMM cells have three
reviewed attempts but only two eligible repeats; their medians remain null.
Full OF/four proteomes remains the sole complete three-eligible-repeat cell.

[Figure v20](threadripper_shared_resource_figure_20261004_v20/shared_threadripper_resources.pdf)
renders the same 19 measured points, one gray-cross exclusion, no imputed
abort point, and only the complete OF cell's median/range. Actual PNG inspected;
all **nine** new reporting tests pass in **1.16s**, including actual PDF bounds,
three-panel colored pixels, source/output hashes, unchanged measured prefix,
missing summaries and refusal of rehabilitated/imputed failures. No accuracy,
scientific setting or isolated efficiency claim changes.

## Resume Boundary

No benchmark job is live and no index-21 request/session exists. Commit/push
the validated checkpoint. New continuation, terminal-review and live-snapshot
sources must use this refreshed lookup and both new resolved sessions. Prepare
prospective current-source recipe/resource/environment/readiness bindings for
index 21, retain original endpoint definitions and shared-host safety bounds,
then check live capacity and actual prepared/native handoff before inference.
Old continuation still points to superseded bindings: do not use it or repeat
its completed preparation. Final panel, manuscript/resource reconciliation
and versioned/archive release remain active and incomplete.
