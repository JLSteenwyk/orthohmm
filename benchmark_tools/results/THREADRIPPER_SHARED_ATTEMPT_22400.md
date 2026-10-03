# Eight-Proteome OrthoFinder Full Run Reviewed

Frozen index 4, full OrthoFinder 3.1.5 on eight proteomes, repeat 0,
completed as job 22400 on the local Threadripper. Parent, batch and native
step all report COMPLETED 0:0. The native command includes sequence search,
MSA/gene-tree inference, species-tree rooting and orthology reconciliation;
this is not the sequence-only checkpoint strategy.

The [retained outcome](threadripper_shared_attempt_22400.json) is an unchanged
copy of the canonical review summary: 28,880 bytes, SHA-256
`eb3a1672f9a11f85d373aa56484ecfaa037070a7257e2f3d3aed4e0ce4646c63`.
The [independent reviewer](review_shared_threadripper_repaired_20261003.py)
passes all four runtime, environment, resource and native-output categories.
Frozen argv/input identity and pre/post runtime checks match; raw resource
counters and whole-run process/pressure observations replay. Underlying raw
evidence stays at the retained work paths referenced by the receipt.

| Endpoint | Retained Observation |
| --- | ---: |
| Input proteins | 165,168 |
| Validated pre-phylogeny checkpoint groups | 33,013 |
| Expanded native ortholog-pair rows | 597,451 |
| Native wall seconds | 1,139.939834238 |
| Native CPU seconds | 17,734.408768 |
| Native-step lifetime peak bytes | 10,987,712,512 |
| Preflight foreign CPU-core equivalents | 52.1288597846 |
| Maximum sampled foreign CPU-core equivalents | 53.6239937821 |

The checkpoint group count validates the input partition; it is not a final
HOG count or accuracy endpoint. Native relation expansion/validation occurs
outside the inference timer. Wall time is the native monotonic command
interval, not the internal log's duration or enclosing scheduler wall. CPU
includes the native-task wrapper bracket; peak is the native-step lifetime
kernel peak including launcher, not process RSS or whole-job teardown memory.
No overhead is subtracted. These are shared-host observations, with unknown,
potentially method-dependent contention distortion, not isolated speedups.

The [five-attempt table](threadripper_shared_panel_snapshot_20261003_v5/panel.json)
and [updated partial figure](threadripper_shared_resource_figure_20261003_v4/shared_threadripper_resources.pdf)
retain four eligible observations, excluded index 0 and all missing repeats.
No cell has three eligible repeats; medians/ranges remain absent. Earlier
snapshots, portable reporting archive and manuscript PDF stay dated and unchanged.

Only after this canonical terminal review does the retained launcher release
index 5, OrthoHMM high sensitivity on eight proteomes, repeat 0, as job 22401.
That submission is not a reviewed result; query its live state on resumption.
Review its terminal evidence before index 6 (full OrthoFinder, 12 proteomes).
Full-panel reporting, manuscript reconciliation and versioned release remain
unfinished. No scientific tuning, retry, quiet-window/DGX requirement or
unrelated job/service change is introduced.
