# Private Threadripper Executor Route

The executor now explicitly selects the retained private deployment rather than
requiring the older shared-runtime plan. This is prospective wiring, not a
production-readiness review, job submission or controlled timing result. Timing
remains deferred; do not stop unrelated work or use the DGX.

## Deployment Contract

New private requests use `deployment: private_v2_20260928` and
`benchmark_tools/run_threadripper_private_scaling.sh` as their exact scheduler
entry point. The private route binds these unchanged artifacts:

| Artifact | SHA256 |
| --- | --- |
| `threadripper_private_commands_20260928.json` | `a9358ac4f3c2f3eb9c1d7ce6a32528f8dd4bffef25d43cf95ac9765c5f2f3d05` |
| `threadripper_private_lookup_v2_20260928.json` | `5996f36dad39c7a7e38c38f134cbe5e13796c4a23dae0444a77a473f12cfdb54` |

The private submission script invokes the already-retained private controller
at `benchmarks/work/threadripper_private_controller_20260928/venv/bin/python`,
after checking its binary SHA256
`8b1cd756be711ef53f35cb6c954472fdfc52094c4619a01f553f64354587388b`.
It has no shared-Python fallback. The executor checks that the bootstrap route,
request route and bound controller entrypoint agree; binary identity alone
does not prove package/library/runtime closure.

Plan and lookup identities propagate through selection, readiness, environment
policy, preflight response, release guard and runtime checker. Shared/private
records cannot be mixed. The recipe must include every top-level benchmark
Python helper and this selected submission script; old recipes are not new
approval. Controller selection does not rewrite runtime manifests or bypass
their source-tree and repeated native-lookup checks.

Requests without a deployment field retain the historical `shared_v3_20260928`
route and original pins/script. That compatibility preserves retained history,
not present shared-runtime validity. Private requests must be explicit; neither
route automatically discovers or upgrades an interpreter or dependency.

## Explicit Current Lookup

The [latest review-helper refresh](results/THREADRIPPER_REVIEW_RUNTIME_REFRESH_20261002.md)
preserves the frozen baseline/controller/plan/private manifest while binding
the five changed and one new Python helpers. All 57,959 records match
before/after native startup, and both native lookups match retained reports.
Actual controller startup covers all 426 observed files and 83 project-module
origins, including the new review combiner. This is startup evidence, not readiness.
The [earlier refresh](results/THREADRIPPER_RUNTIME_REFRESH_20261002.md) and its
327-case test receipt remain historical evidence, not repeated tests.

An explicitly private request may supply optional `runtime_lookup`, a direct
absolute file record with exactly `path`, positive integer `bytes` and lowercase
64-digit `sha256`. Its digest must equal `lookup_sha256` in both request and
readiness review. Shared/implicit routes reject it. The old default is untouched.
This avoids embedding a refreshed lookup hash into source inventoried by that
lookup. No runtime discovery, fallback or implicit readiness approval is added.

The [latest local-host lookup receipt](results/threadripper_private_lookup_review_20261002.json)
has 6,584 bytes and SHA256
`be7084d73bd90dedcdec210bbc3f9575eb9e3d9461c5f882d7b958c0c3b670de`.
Bind its actual direct local path in `runtime_lookup`; it refers to retained
Threadripper work artifacts, not a portable runtime distribution. Do not
calculate replacement hashes from untrusted manifests or reuse old readiness.
All other request fields, full current-source recipe, history, real policy,
observer validation and environmental handoff remain required.

The [latest actual controller-startup check](results/THREADRIPPER_REVIEW_RUNTIME_REFRESH_20261002.md)
also imports the current process/pressure and panel reviewers:
412 file-backed modules, 83 project modules and 426 observed files all match
pinned identities. Pure production/engineering deployment resolution passes. This is not
executor selection, real environmental handoff or complete workload/file-I/O
closure; no scheduler job or measurement function executes.

## Preserved Gates

Private execution now requires `threadripper_environment_policy_v2` with
explicit `native_pressure_role: diagnostic_only`. Whole-native-interval PSI
magnitudes are retained diagnostics, not exclusion criteria for a method's
own stalls. Evidence integrity, strict outside-process/CPU policy and parked
prelaunch pressure limits still gate execution. Historical shared/v1 replay
semantics remain unchanged. See the [prospective correction and retained
validation](results/THREADRIPPER_PRESSURE_ROLES_20261002.md). No real approved
policy or live handoff is created; old source/runtime recipes are not new
readiness approval after these helper changes.

All 27 production identities, datasets, scientific settings, order, input
preparation, 32-CPU native affinity, 64 scheduler slots, 128 GiB allocation,
26-hour limit and no-requeue policy remain unchanged. The private script does
not call `sbatch`, release held jobs, retry or authorize subsequent runs.

The executor still requires source/readiness/history reviews, fresh same-boot
environment handoff, live allocation budget, runtime checks, whole-run process
and pressure review, and independent post-run native/resource/environment
admission. A synthetic test's review fields are not evidence of those real
requirements. A failed controller check precedes session creation and native
work. Every actual submission/failed bootstrap must remain in the external
attempt ledger; it is not silently retried.

## Remaining Work

The [conditional panel-review workflow](results/THREADRIPPER_OVERHEAD_PANEL_REVIEW_20261002.md)
now combines native replay with all four recorded engineering review categories.
It distinguishes complete reviewed numerical pass/failure from missing evidence
and never authorizes production or independently certifies review conclusions.
All 137 focused cases pass. The actual CLI preserves all 54 tasks as unrun
with a null budget; no real passing review or timing result is created.

The [boundary post-run correction](results/THREADRIPPER_BOUNDARY_ENVIRONMENT_REVIEW_20261002.md)
passes the validated task arm to environmental review. Boundary-only PSI now
has an explicit two-point diagnostic scope; it is not required to satisfy
periodic native-point cadence. Whole-run process monitoring and all default
periodic/production behavior remain unchanged. All 643 focused cases pass;
offline synthetic decisions are fail/pass/fail, not real host eligibility.

The [native-observer control route](results/THREADRIPPER_OVERHEAD_EXECUTOR_20261002.md)
now selects either collector for the separate 54-task engineering panel.
It requires an explicit private lookup and `overhead_plan`, engineering
readiness/history schemas and the same environmental safeguards. Failure
stops continuation; no engineering result becomes a production identity.
These helper bytes are covered by the latest review-helper inventory/startup
receipt above, not by the earlier inventory. Recheck the bound runtime/source
at launch without changing the scientific plan or repeating calibration.
Actual overhead/environmental execution is pending.

The passed 22380 calibration remains accounting/cadence evidence, not causal
slowdown. The separate 54-task native observer-control panel and 27 production
identities remain unadmitted. A real reviewed process/service policy, complete
native environmental handoff, launch-time current-source/runtime validation and
quiet window are still required. A retained lookup pin is not proof that its
historical transitive runtime files still match after later helper changes.

This route has synthetic composition/negative checks, Bash syntax validation
and the bounded actual private-controller inventory/native-import check above,
not full executor/environmental integration. The original submission script,
scientific package, historical manifests and diagnostic receipts are unchanged.
See the [milestone evidence](results/THREADRIPPER_PRIVATE_EXECUTOR_20261002.md).
