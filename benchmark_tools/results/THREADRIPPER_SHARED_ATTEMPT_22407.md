# Four-Proteome High-Sensitivity Repeat Reviewed

Frozen index 11, OrthoHMM high sensitivity on four proteomes, repeat 1,
completed as job 22407 on the local Threadripper. Parent, batch and native
Slurm step report COMPLETED 0:0. Parent elapsed is 9:14; native step elapsed
is 6:43. HMM search, threshold/graph construction and clustering finish under
the frozen configuration, without phylogenetic reconciliation. Native exit
alone was not treated as enclosing-job completion.

The [retained outcome](threadripper_shared_attempt_22407.json) is an unchanged
copy of the canonical review summary: 4,466 bytes, SHA-256
`8a13bd4baa02e2ce51e87010d8a0a2974960dd6195d9b1bbe50dcd45f70fda40`.
The [independent reviewer](review_shared_threadripper_repaired_20261003.py)
passes runtime, environment, resources and native-output checks. Readback
verifies all four category pins and decisions. Frozen argv/input identities,
pre/post runtime bindings, raw resource counters and whole-run sampled
process/pressure observations agree with the retained evidence.

| Endpoint | Retained Observation |
| --- | ---: |
| Input proteins | 73,266 |
| Output orthogroups | 35,242 |
| Native wall seconds | 377.661928371 |
| Native CPU seconds | 9,377.291271 |
| Native-step lifetime peak bytes | 3,387,318,272 |
| Preflight foreign CPU-core equivalents | 54.4996349095 |
| Maximum sampled foreign CPU-core equivalents | 54.5505774455 |

Output checks validate a complete native partition, not prediction accuracy.
Wall time is the independently replayed native monotonic command interval.
CPU includes the native-task-subtree wrapper bracket; memory is native-step
lifetime peak including launcher, not process RSS or whole-job teardown memory.
Process and pressure failure maps are empty. Contention distortion remains
unknown and potentially method dependent; no overhead subtraction or isolated
ranking is claimed. Preserve this repeat without fastest-repeat selection.

The [twelve-attempt table](threadripper_shared_panel_snapshot_20261003_v12/panel.json)
and [updated partial figure](threadripper_shared_resource_figure_20261003_v11/shared_threadripper_resources.pdf)
retain eleven eligible observations and excluded index 0. This is a scheduled
repeat, not a retry of that failure. The four-proteome high-sensitivity cell
has one eligible repeat and one excluded attempt. No cell has three eligible
repeats; all medians/ranges remain unavailable. Earlier snapshots, reporting
archive and manuscript PDF retain their unchanged dated payloads.

Only after canonical terminal review does the launcher release index 12,
full OrthoFinder 3.1.5 on eight proteomes, repeat 1, as job 22408. The
[live checkpoint](threadripper_shared_live_22408.json) is dated evidence, not a
terminal review; query that exact job on resumption. Review its terminal evidence
before index 13 (OrthoHMM high sensitivity, eight proteomes, repeat 1).
Full-panel reporting, manuscript reconciliation and versioned release remain
unfinished. No scientific tuning, retry, quiet-window/DGX requirement or
unrelated workload change is introduced.
