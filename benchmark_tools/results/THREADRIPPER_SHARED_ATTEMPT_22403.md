# Twelve-Proteome OrthoHMM High-Sensitivity Run Reviewed

Frozen index 7, OrthoHMM high sensitivity on 12 proteomes, repeat 0,
completed as job 22403 on the local Threadripper. Parent, batch and native
step all report COMPLETED 0:0. Built-in HMM search, threshold/graph construction
and clustering complete with the unchanged frozen method. Phylogenetic
reconciliation is not enabled in this configuration.

The [retained outcome](threadripper_shared_attempt_22403.json) is an unchanged
copy of the canonical review summary: 4,471 bytes, SHA-256
`beb31a9d114d41a18f2e8df81f82ec637ee9a7d46269a882a5c5745da127ef23`.
The [independent reviewer](review_shared_threadripper_repaired_20261003.py)
passes runtime, environment, resources and native-output checks. Frozen
argv/input identity and pre/post runtime checks match; raw resource counters
and whole-run process/pressure observations independently replay. Underlying
evidence remains at retained local paths in the receipt.

| Endpoint | Retained Observation |
| --- | ---: |
| Input proteins | 251,378 |
| Output orthogroups | 62,885 |
| Native wall seconds | 3,027.060520509 |
| Native CPU seconds | 86,415.319181 |
| Native-step lifetime peak bytes | 11,309,355,008 |
| Preflight foreign CPU-core equivalents | 58.8194471963 |
| Maximum sampled foreign CPU-core equivalents | 60.3531642451 |

Output checks validate a complete native partition, not prediction accuracy.
Wall time is the native monotonic command interval. CPU includes the native
task-subtree wrapper bracket; memory is the native-step lifetime kernel peak
including launcher, not process RSS or whole-job teardown memory. No overhead
is subtracted. Shared-host contention may distort timing by an unknown,
potentially method-dependent amount; no isolated speedup or ranking follows.

The [eight-attempt table](threadripper_shared_panel_snapshot_20261003_v8/panel.json)
and [updated partial figure](threadripper_shared_resource_figure_20261003_v7/shared_threadripper_resources.pdf)
retain seven eligible observations and excluded index 0. High sensitivity and
full OrthoFinder each have one eligible 12-proteome observation, but no cell
has three eligible repeats. No median/range is imputed. Earlier snapshots,
reporting archive and manuscript PDF retain their unchanged dated payloads.

Only after canonical terminal review does the launcher release index 8,
OrthoHMM satellite_v2 on 12 proteomes, repeat 0, as job 22404. That submission
is not a reviewed result; query live state on resumption. Review its terminal
evidence before index 9 (satellite_v2, four proteomes, repeat 1).
Full-panel reporting, manuscript reconciliation and versioned release remain
unfinished. No scientific tuning, retry, quiet-window/DGX requirement or
unrelated job/service change is introduced.
