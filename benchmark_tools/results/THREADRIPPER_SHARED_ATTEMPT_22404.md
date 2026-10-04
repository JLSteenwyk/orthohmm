# Twelve-Proteome OrthoHMM Phylogeny Run Reviewed

Frozen index 8, OrthoHMM satellite_v2 with inferred phylogeny on 12 proteomes,
repeat 0, completed as job 22404 on the local Threadripper. Parent, batch and
native Slurm step report COMPLETED 0:0. Built-in HMM search, graph construction,
grouping, species-tree/gene-tree inference and reconciliation finish with the
unchanged frozen configuration. Native exit alone was not treated as completion.

The [retained outcome](threadripper_shared_attempt_22404.json) is an unchanged
copy of the canonical review summary: 5,167 bytes, SHA-256
`15e34be890d8a86335886bfb902e4e9b882a342da9ae5111231c902ff80c8f10`.
The [independent reviewer](review_shared_threadripper_repaired_20261003.py)
passes runtime, environment, resources and native-output checks. Frozen
argv/input identities and pre/post runtime bindings match; raw resource
counters and whole-run process/pressure observations independently replay.
Underlying evidence remains at the retained local paths in the receipt.

| Endpoint | Retained Observation |
| --- | ---: |
| Input proteins | 251,378 |
| Output orthogroups / root HOGs | 59,770 |
| Native pair rows | 966,439 |
| Native wall seconds | 4,809.292429717 |
| Native CPU seconds | 136,293.332246 |
| Native-step lifetime peak bytes | 11,285,721,088 |
| Preflight foreign CPU-core equivalents | 59.9950600778 |
| Maximum sampled foreign CPU-core equivalents | 61.1453857480 |

Output checks validate a complete native partition and native relations, not
prediction accuracy. Wall time is the native monotonic command interval.
CPU includes the native-task-subtree wrapper bracket; memory is the native-step
lifetime kernel peak including launcher, not process RSS or whole-job teardown
memory. No overhead is subtracted. Contention may distort timing by an unknown,
potentially method-dependent amount; no isolated speedup or ranking follows.

The [nine-attempt table](threadripper_shared_panel_snapshot_20261003_v9/panel.json)
and [updated partial figure](threadripper_shared_resource_figure_20261003_v8/shared_threadripper_resources.pdf)
retain eight eligible observations and excluded index 0. Each method now has
one eligible 12-proteome observation, but no cell has three eligible repeats.
No median/range is imputed. Earlier snapshots, reporting archive and manuscript
PDF retain their unchanged dated payloads.

Only after canonical terminal review does the launcher release index 9,
OrthoHMM satellite_v2 on four proteomes, repeat 1, as job 22405. That submission
is not a reviewed result; query live state on resumption. Review its terminal
evidence before index 10 (full OrthoFinder 3.1.5, four proteomes, repeat 1).
Full-panel reporting, manuscript reconciliation and versioned release remain
unfinished. No scientific tuning, retry, quiet-window/DGX requirement or
unrelated job/service change is introduced.
