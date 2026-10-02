# Private Deployment Executor Preparation

The executor's historical default still selected the shared plan/lookup and
the submission script still invoked shared Anaconda. The retained private
deployment had separately validated fixtures, but no explicit executor route.
This milestone wires that route without updating historical manifests or
creating a passing readiness review, request or scheduler job.

See the [private request and bootstrap contract](../THREADRIPPER_PRIVATE_EXECUTOR.md)
and [machine-readable review](threadripper_private_executor_review_20261002.json).
The review is 7,491 bytes, SHA256
`f9ec6bf4145068824128ffac129fc219d0a88e93bd61ad339d7ef92dc0a52669`.
Base source is `84d4e9bdfc026c10b6cdb170d4467c85daaa3e34`.

Private requests explicitly select `private_v2_20260928`, the unchanged private
27-run plan and v2 lookup pins, and the new private submission script. The
private controller's binary is checked before bootstrap and its literal
entrypoint is checked against the lookup binding before session creation.
No shared-interpreter fallback or arbitrary plan/interpreter discovery exists.
Selection propagates through policy, readiness, responder, release guard and
runtime checker. Omitted deployment fields retain historical shared behavior,
not current shared-runtime approval. Unknown or mixed deployment records fail.

Six retained artifact files and the original submission script are byte-identical
to the base commit. The bound private controller is a direct file, not a symlink;
its actual binary pin matches. This is static binary identity, not current
package/library/transitive runtime closure or a fresh native lookup. Scientific
package, settings, scores, 27 identities, resources/order and no-retry gates are
unchanged. No timing, host contention poll, DGX access or unrelated process/
service action occurs.

## Actual Checks

Initial existing route/worker/lifecycle panel: 81 passes. First expanded panel
has 91 passes/10 failures in 7.15s: two synthetic plans had identical bytes,
and eight private-controller fixtures bound the resolved base-Python path.
After differentiating plan metadata and preserving the virtual-environment
path, the second expanded panel has 114 passes/eight failures in 7.46s: the
strict existing pin checker rejects symlinked synthetic interpreter records.
Both failed XML receipts remain retained. The actual bound private controller
is a direct file and satisfies the unchanged pin checker.

Use a direct, explicitly synthetic controller file for composition tests;
its content is never executed. Do not weaken controller path/digest assertions.
The corrected three-module panel passes 122 cases in 7.04s. Add missing/wrong-
path/wrong-hash controller rejection and broaden unchanged history, progress,
stream and runtime-checker regression. Final **229 cases pass in 8.34s**, zero
failures/errors/skips. JUnit is 34,951 bytes, SHA256
`657788efa8f3b05e4521e34677b395c2c36fdfd571e9347b5714b95b328a6a98`.
Panels overlap, not additive. Bash syntax validation passes; the submission
script is not executed. Native measurement and policy conclusions are mocked
or explicitly synthetic in these composition tests, not real readiness evidence.

At 14:27:25 UTC exact prior-source CI run 37017458701 has successful wheel/Linux
jobs, live five macOS jobs and queued docs. No execution counts are downloaded
or inferred; this run excludes the new private-route source. Post-push CI must
be observed separately.

## Remaining Boundary

No production permit or final recipe/readiness artifact is generated. The
new route still needs final current-source/runtime verification, a real reviewed
environment policy, actual private-controller/native handoff validation, causal
observer slowdown assessment and a quiet window. The passed 22380 accounting/
cadence calibration is reused, not rerun or upgraded to slowdown evidence.
Historical runtime bindings may still reject later helper drift; they are not
silently refreshed by this change. The 54 engineering tasks are separate from
the 27 production identities, and no attempt is launched here.

Controlled resource evidence, remaining QfO uncertainty, complete release/
rights/deposition and final manuscript/figure/archive reconciliation remain
unfinished. The publication goal remains active.
