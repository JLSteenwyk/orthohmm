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

## Preserved Gates

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

The passed 22380 calibration remains accounting/cadence evidence, not causal
slowdown. The separate 54-task native observer-control panel and 27 production
identities remain unadmitted. A real reviewed process/service policy, complete
native environmental handoff, final current-source/runtime validation and
quiet window are still required. A retained lookup pin is not proof that its
historical transitive runtime files still match after later helper changes.

This route has synthetic composition/negative checks and Bash syntax validation,
not a new native/private-controller integration run. The original submission
script, scientific package, manifests and diagnostic receipts are unchanged.
See the [milestone evidence](results/THREADRIPPER_PRIVATE_EXECUTOR_20261002.md).
