# Twelve-Proteome High-Sensitivity Repeat Reviewed

Frozen index 15, OrthoHMM high sensitivity on twelve proteomes, repeat 1,
completed as job 22411 on the local Threadripper. Parent, batch and native
Slurm step report COMPLETED 0:0. Parent elapsed is 49:45; native step elapsed
is 45:48. Built-in HMM search, threshold estimation, graph construction and
clustering finish under the frozen configuration. This configuration stops
after grouping, without phylogenetic reconciliation. Worker/main CPU counters
continued increasing during unchanged progress displays. Exited search or
clustering workers were not treated as stopped enclosing inference. Native
exit alone was not treated as enclosing-job completion.

The [retained outcome](threadripper_shared_attempt_22411.json) is an unchanged
copy of the canonical review summary: 4,473 bytes, SHA-256
`14fb0e58a31d2ac4278871383e7b56ae78029c544a1121d4673507247626513f`.
The [independent reviewer](review_shared_threadripper_repaired_20261003.py)
passes runtime, environment, resources and native-output checks. Readback
verifies all four category pins and decisions. Frozen argv/input identities,
pre/post runtime bindings, raw resource counters and whole-run sampled
process/pressure observations agree with the retained evidence.

| Endpoint | Retained Observation |
| --- | ---: |
| Input proteins | 251,378 |
| Orthogroups, including singletons | 62,885 |
| Native wall seconds | 2,720.077391309 |
| Native CPU seconds | 78,335.830623 |
| Native-step lifetime peak bytes | 11,521,089,536 |
| Preflight foreign CPU-core equivalents | 54.0887580242 |
| Maximum sampled foreign CPU-core equivalents | 54.0032285143 |

The complete partition is not a new accuracy evaluation. Wall time is the
independently replayed native monotonic command interval. CPU includes the
native-task-subtree wrapper bracket; memory is native-step lifetime peak
including launcher, not process RSS or whole-job teardown memory. Process and
pressure failure maps are empty. Preflight and whole-run sampled foreign
demand use different observation windows; their values are not a slowdown
correction. Contention distortion remains unknown and potentially method
dependent; no overhead subtraction or isolated ranking is claimed. Preserve
this repeat without fastest-repeat selection.

The [sixteen-attempt table](threadripper_shared_panel_snapshot_20261003_v16/panel.json)
and [updated partial figure](threadripper_shared_resource_figure_20261003_v15/shared_threadripper_resources.pdf)
retain fifteen eligible observations and excluded index 0. High sensitivity's
twelve-proteome cell has two eligible repeats; no cell has three eligible
repeats and all medians/ranges remain unavailable. Earlier snapshots, reporting
archive and manuscript PDF retain their unchanged dated payloads. The
[eight-proteome phylogenetic repeat difference](THREADRIPPER_SHARED_ATTEMPT_22410.md)
remains documented separately; it is not erased by this successful run.

Only after canonical terminal review does the launcher release index 16,
OrthoHMM satellite_v2 on twelve proteomes, repeat 1, as job 22412. The
[live checkpoint](threadripper_shared_live_22412.json) is dated evidence, not a
terminal review; query that exact job on resumption. Review its terminal
evidence before index 17 (full OrthoFinder 3.1.5, twelve proteomes, repeat 1).
Full-panel reporting, manuscript reconciliation and versioned release remain
unfinished. No scientific tuning, retry, quiet-window/DGX requirement or
unrelated workload change is introduced.
