# Eight-Proteome Full OrthoFinder Repeat Reviewed

Frozen index 12, full OrthoFinder 3.1.5 on eight proteomes, repeat 1, completed
as job 22408 on the local Threadripper. Parent, batch and native Slurm step
report COMPLETED 0:0. Parent elapsed is 23:06; native step elapsed is 20:07.
All 64 DIAMOND searches, MCL grouping, 1,828 species-tree alignments, 12,110
remaining alignment/tree tasks, STRIDE rooting, reconciliation and final writing
finish under the frozen configuration. The native log retains MCL's warning
about 22 overlap instances. Long FastTree/FAMSA processes accumulated CPU;
stable file counts were not treated as a stopped job. Native exit alone was
not treated as enclosing-job completion.

The [retained outcome](threadripper_shared_attempt_22408.json) is an unchanged
copy of the canonical review summary: 28,882 bytes, SHA-256
`4dec1a8b5d37d4612c175fb32ec84b14d1d525890150d22f71124b0d18745bc7`.
The [independent reviewer](review_shared_threadripper_repaired_20261003.py)
passes runtime, environment, resources and native-output checks. Readback
verifies all four category pins and decisions. Frozen argv/input identities,
pre/post runtime bindings, raw resource counters and whole-run sampled
process/pressure observations agree with the retained evidence.

| Endpoint | Retained Observation |
| --- | ---: |
| Input proteins | 165,168 |
| Pre-phylogeny checkpoint groups, including singletons | 33,013 |
| Expanded native ortholog-pair rows | 597,451 |
| Native wall seconds | 1,179.970330278 |
| Native CPU seconds | 18,787.311821 |
| Native-step lifetime peak bytes | 10,967,740,416 |
| Preflight foreign CPU-core equivalents | 54.2986048129 |
| Maximum sampled foreign CPU-core equivalents | 55.0811204187 |

Checkpoint groups are not a phylogenetic pair prediction or an accuracy score.
Native pair validation expands relation rows outside the inference timer.
Wall time is the independently replayed native monotonic command interval,
not the 1,174.491257 seconds printed in OrthoFinder's log. CPU includes the
native-task-subtree wrapper bracket; memory is native-step lifetime peak
including launcher, not process RSS or whole-job teardown memory. Process and
pressure failure maps are empty. Contention distortion remains unknown and
potentially method dependent; no overhead subtraction or isolated ranking is
claimed. Preserve this repeat without fastest-repeat selection.

The [thirteen-attempt table](threadripper_shared_panel_snapshot_20261003_v13/panel.json)
and [updated partial figure](threadripper_shared_resource_figure_20261003_v12/shared_threadripper_resources.pdf)
retain twelve eligible observations and excluded index 0. Full OrthoFinder's
eight-proteome cell has two eligible repeats; no cell has three eligible
repeats and all medians/ranges remain unavailable. Earlier snapshots, reporting
archive and manuscript PDF retain their unchanged dated payloads.

Only after canonical terminal review does the launcher release index 13,
OrthoHMM high sensitivity on eight proteomes, repeat 1, as job 22409. The
[live checkpoint](threadripper_shared_live_22409.json) is dated evidence, not a
terminal review; query that exact job on resumption. Review its terminal
evidence before index 14 (OrthoHMM satellite_v2, eight proteomes, repeat 1).
Full-panel reporting, manuscript reconciliation and versioned release remain
unfinished. No scientific tuning, retry, quiet-window/DGX requirement or
unrelated workload change is introduced.
