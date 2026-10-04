# Phylogenetic Repeat Aborted Before Native Inference

After index 19 review and reporting commit/push `7de68449`, the unchanged
continuation helper validates the twenty-identity prefix and selects/holds/
request-binds/releases only index **20** as **22416**: phylogenetic OrthoHMM
satellite_v2, four proteomes, repeat 2, 73,266 frozen input proteins. The
[public submission](threadripper_shared_submission_22416.json) matches canonical
launch bytes: 1,175 bytes, SHA-256
`94726ed6aff439f459692e2531f4e1b9360bfb99065e347a1d51b233674cc254`.
Request is 7,157 bytes, SHA-256
`c72c9ab857a8a70a72da1fe8b6dd18b2c7074b0ff43b579990b62311ce61faf0`.
Frozen science/input/runtime/readiness and 32-physical-CPU/128-GiB limits are
unchanged. No excluded identity is retried or unrelated workload displaced.

Fresh accounting reports parent/batch FAILED 1:0 at 3:06. The parked native
step completed 0:0 at 0:38, but the gate received `{"abort": true}` and wrote
`observer_did_not_release_native`. OrthoHMM never started. There is no native
log, completion record, pressure-point series, metrics file or nonempty output.
The wrapper prepared an empty output directory; that is not native output.
Native wall/CPU/peak-memory endpoints are absent, not zero or imputed.

## Retained Failure Review

The [failure audit](threadripper_shared_prenative_failure_22416.json), 8,526
bytes, SHA-256
`53bac27059a4a636b3986d98ad57907fb106899b645cf103ce6db1f9b3bf507e`,
binds request/plan, actual abort gate, preflight detail, worker lifecycle,
pre/post recorded runtime checks, immutable first sample, final two-sample
stream, controller and fresh accounting. The
[terminal controller](threadripper_terminal_controller_22416.json) is captured
before purge with actual timestamps and matches its canonical copy: 1,830
bytes, SHA-256
`41bb9a9c05d88747461a503f6aec85b47c79514f0960fa10a21b3d4c5c8fc43f`.
The audit does not turn the failed preflight into a passing resource review.

The wrapper records `Environmental worker response deadline exceeded`.
Observer preparation is not synchronized with request publication: source
inspection shows deployment/history/policy/source checks before waiting for
the marker. Actual review starts 11.799922912s after request; its response is
published 21.377395496s after request, beyond the frozen 20s window.
Final observer detail records missing `/proc/1512928/cgroup`, identifying the
parked worker, and `deadline_expired`. The wrapper had already taken its
timeout/abort path. These records support a deadline/handoff failure, not a
diagnosis of OrthoHMM performance or a quantitative explanation of slow setup.
Context cleanup subsequently reports the worker was not successfully joined;
both errors remain preserved.

Available memory is 329,335,861,248 / 329,492,717,568 bytes, both above the
137,438,953,472-byte floor. Background-process policy review passes and records
40.3849069388 foreign CPU-core equivalents. Those diagnostics do not establish
a passed preflight, complete workload attribution, isolation or causal slowdown.
The immutable-sample repair remains in place; this is not the earlier mutable
stream hash failure.

The first audit implementation incorrectly rejected the empty preparation
directory and produced no receipt. Correct only artifact classification, not
the failed production outcome. Twenty initial assessment tests pass in 3.88s;
after adding actual inventory controls, all 27 focused tests pass in 3.79s.
They reject contradictory identities/deadlines/gates, unsafe capacity, changed
initial observations, any native artifacts and output symlinks. Successful
audit then independently checks retained evidence once; no native rerun.

## Remaining Work

No benchmark job is live. No index-20 successful canonical summary, live
inference snapshot, retry or index-21 submission is fabricated. Existing
table v20/figure v19 remain the twenty-identity resolved-prefix checkpoint,
predating this separately audited failure; their index-20 pending row is not
a statement that the attempt remains unsubmitted.

Repair preparation/request synchronization prospectively while preserving
the fresh-capacity check, 20s release response and 45s parked-worker gate,
monitoring cadence, immutable evidence, runtime integrity and native budget.
Validate the repair and refresh affected source/runtime/protocol bindings;
explicitly resolve this excluded history before index 21. Extend reporting
and continuation to retain both pre-native aborts without inventing native
endpoints or retrying either. Do not bypass current continuation's required
successful index-20 summary. Final panel, manuscript/resource reconciliation
and versioned/archive release remain active and incomplete.
