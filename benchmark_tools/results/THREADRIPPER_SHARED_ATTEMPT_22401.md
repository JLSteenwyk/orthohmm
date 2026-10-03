# Eight-Proteome OrthoHMM High-Sensitivity Run Reviewed

Frozen index 5, OrthoHMM high sensitivity on eight proteomes, repeat 0,
completed as job 22401 on the local Threadripper. Parent, batch and native
step all report COMPLETED 0:0. Search, threshold/graph construction and
clustering complete using the unchanged frozen method; phylogenetic
reconciliation is not enabled in this configuration.

The [retained outcome](threadripper_shared_attempt_22401.json) is an unchanged
copy of the canonical review summary: 4,469 bytes, SHA-256
`a500f67ca948ad4cef14b3903e09b301ce0ef892973698ce7d306666e909544e`.
The [independent reviewer](review_shared_threadripper_repaired_20261003.py)
passes runtime, environment, resources and native-output checks. Frozen
argv/input identity and pre/post runtime checks match; raw resource counters
and whole-run process/pressure observations independently replay. Underlying
evidence remains at retained local paths in the receipt.

| Endpoint | Retained Observation |
| --- | ---: |
| Input proteins | 165,168 |
| Output orthogroups | 58,278 |
| Native wall seconds | 1,238.935206005 |
| Native CPU seconds | 34,395.748088 |
| Native-step lifetime peak bytes | 8,293,474,304 |
| Preflight foreign CPU-core equivalents | 52.0202008983 |
| Maximum sampled foreign CPU-core equivalents | 52.9808092544 |

Output checks validate a complete native partition, not prediction accuracy.
Wall time is the native monotonic command interval. CPU includes the native
task-subtree wrapper bracket; memory is the native-step lifetime kernel peak
including launcher, not process RSS or whole-job teardown memory. No overhead
is subtracted. Shared-host contention may distort timings by an unknown,
potentially method-dependent amount; no isolated speedup or ranking follows.

The [six-attempt table](threadripper_shared_panel_snapshot_20261003_v6/panel.json)
and [updated partial figure](threadripper_shared_resource_figure_20261003_v5/shared_threadripper_resources.pdf)
retain five eligible observations and excluded index 0. Each method now has
one eligible eight-proteome observation, but no cell has three eligible
repeats. No median/range is imputed. Earlier snapshots, reporting archive and
manuscript PDF remain unchanged dated checkpoints.

Only after the canonical terminal review does the launcher release index 6,
full OrthoFinder 3.1.5 on 12 proteomes, repeat 0, as job 22402. That submission
is not a reviewed result; check live state on resumption. Review its terminal
evidence before index 7 (OrthoHMM high sensitivity on 12 proteomes).
Full-panel reporting, manuscript reconciliation and versioned release remain
unfinished. No scientific tuning, retry, quiet-window/DGX requirement or
unrelated job/service change is introduced.
