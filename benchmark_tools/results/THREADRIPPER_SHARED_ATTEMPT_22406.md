# Four-Proteome Full OrthoFinder Repeat Reviewed

Frozen index 10, full OrthoFinder 3.1.5 on four proteomes, repeat 1, completed
as job 22406 on the local Threadripper. Parent, batch and native Slurm step
report COMPLETED 0:0. Parent elapsed is 11:11; native step elapsed is 8:20.
All 16 DIAMOND searches, MCL grouping, 1,771 species-tree alignments, 10,466
remaining alignment/tree tasks, STRIDE rooting, reconciliation and final writing
finish under the frozen configuration. Native exit alone was not treated as
enclosing-job completion.

The [retained outcome](threadripper_shared_attempt_22406.json) is an unchanged
copy of the canonical review summary: 10,329 bytes, SHA-256
`628b0a3fd050a0ade748be71f768c9a462c8cac728200a615b8ecc08c4ab411d`.
The [independent reviewer](review_shared_threadripper_repaired_20261003.py)
passes runtime, environment, resources and native-output checks. Readback
verifies all four category pins and decisions. Frozen argv/input identities,
pre/post runtime bindings, raw resource counters and whole-run sampled
process/pressure observations agree with the retained evidence.

| Endpoint | Retained Observation |
| --- | ---: |
| Input proteins | 73,266 |
| Pre-phylogeny checkpoint groups, including singletons | 24,052 |
| Expanded native ortholog-pair rows | 88,890 |
| Native wall seconds | 474.299930051 |
| Native CPU seconds | 5,008.896961 |
| Native-step lifetime peak bytes | 7,268,642,816 |
| Preflight foreign CPU-core equivalents | 53.0266716099 |
| Maximum sampled foreign CPU-core equivalents | 53.9775442496 |

Checkpoint groups are not a phylogenetic pair prediction or an accuracy score.
Native pair validation expands relation rows outside the inference timer.
Wall time is the independently replayed native monotonic command interval,
not the 468.843407 seconds printed in OrthoFinder's log. CPU includes the
native-task-subtree wrapper bracket; memory is native-step lifetime peak
including launcher, not process RSS or whole-job teardown memory. Process and
pressure failure maps are empty. Contention distortion remains unknown and
potentially method dependent; no overhead subtraction or isolated ranking is
claimed. Preserve this repeat without fastest-repeat selection.

The [eleven-attempt table](threadripper_shared_panel_snapshot_20261003_v11/panel.json)
and [updated partial figure](threadripper_shared_resource_figure_20261003_v10/shared_threadripper_resources.pdf)
retain ten eligible observations and excluded index 0. Both four-proteome
phylogenetic cells have two eligible repeats, not complete three-repeat
summaries. All cell medians/ranges remain unavailable. Earlier snapshots,
reporting archive and manuscript PDF retain their unchanged dated payloads.

Only after canonical terminal review does the launcher release index 11,
OrthoHMM high sensitivity on four proteomes, repeat 1, as job 22407. The
[live checkpoint](threadripper_shared_live_22407.json) is dated evidence, not a
terminal review; query that exact job on resumption. Review its terminal
evidence before index 12 (full OrthoFinder 3.1.5, eight proteomes, repeat 1).
Full-panel reporting, manuscript reconciliation and versioned release remain
unfinished. No scientific tuning, retry, quiet-window/DGX requirement or
unrelated workload change is introduced.
