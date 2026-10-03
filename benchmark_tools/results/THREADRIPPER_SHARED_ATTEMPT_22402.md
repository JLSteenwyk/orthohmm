# Twelve-Proteome Full OrthoFinder Run Reviewed

Frozen index 6, full OrthoFinder 3.1.5 on 12 proteomes, repeat 0, completed
as job 22402 on the local Threadripper. Parent, batch and native step all
report COMPLETED 0:0. The native log confirms all 144 DIAMOND searches,
1,364 species-tree alignments and 13,342 remaining alignment/tree tasks,
followed by STRIDE, reconciliation and hierarchical-group output. This is
the full phylogenetic pipeline, not a sequence-only checkpoint run.

The [retained outcome](threadripper_shared_attempt_22402.json) is an unchanged
copy of the canonical review summary: 60,845 bytes, SHA-256
`ca6d16ea6e1ebb58b9c8c3779f5240f8381adcbb5f9ce00ca075ec95959c8b7c`.
The [independent reviewer](review_shared_threadripper_repaired_20261003.py)
passes runtime, environment, resources and native-output checks. Frozen
argv/input identity and pre/post runtime checks match; raw resource counters
and whole-run process/pressure observations independently replay. Underlying
evidence remains at retained local paths in the receipt.

| Endpoint | Retained Observation |
| --- | ---: |
| Input proteins | 251,378 |
| Pre-phylogeny checkpoint groups | 34,230 |
| Expanded native ortholog pair rows | 1,487,084 |
| Native wall seconds | 1,896.732060134 |
| Native CPU seconds | 39,246.594856 |
| Native-step lifetime peak bytes | 11,897,724,928 |
| Preflight foreign CPU-core equivalents | 52.7452027652 |
| Maximum sampled foreign CPU-core equivalents | 60.8152990336 |

Checkpoint groups are not final hierarchical orthogroups. Native relation-row
expansion is output validation outside the inference timer, not a new accuracy
score. Wall time is the native monotonic command interval. CPU includes the
native task-subtree wrapper bracket; memory is the native-step lifetime kernel
peak including launcher, not process RSS or whole-job teardown memory. No
overhead is subtracted. Shared-host contention may distort timing by an unknown,
potentially method-dependent amount; no isolated speedup or ranking follows.

The [seven-attempt table](threadripper_shared_panel_snapshot_20261003_v7/panel.json)
and [updated partial figure](threadripper_shared_resource_figure_20261003_v6/shared_threadripper_resources.pdf)
retain six eligible observations and excluded index 0. OrthoFinder has one
eligible observation at each size, not three-repeat summaries. No median/range
is imputed. Earlier snapshots, reporting archive and manuscript PDF retain
their unchanged dated payloads.

Only after canonical terminal review does the launcher release index 7,
OrthoHMM high sensitivity on 12 proteomes, repeat 0, as job 22403. That
submission is not a reviewed result; query live state on resumption. Review
its terminal evidence before index 8 (OrthoHMM satellite_v2, 12 proteomes).
Full-panel reporting, manuscript reconciliation and versioned release remain
unfinished. No scientific tuning, retry, quiet-window/DGX requirement or
unrelated job/service change is introduced.
