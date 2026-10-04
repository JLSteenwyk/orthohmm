# Eight-Proteome Phylogenetic Repeat Reviewed

Frozen index 14, OrthoHMM satellite_v2 on eight proteomes, repeat 1,
completed as job 22410 on the local Threadripper. Parent, batch and native
Slurm step report COMPLETED 0:0. Parent elapsed is 39:03; native step elapsed
is 35:15. Built-in HMM search, threshold estimation, graph construction,
clustering, species-tree inference, family alignment/tree inference and
reconciliation finish under the frozen configuration. Native exit alone
was not treated as enclosing-job completion.

The [retained outcome](threadripper_shared_attempt_22410.json) is an unchanged
copy of the canonical review summary: 5,165 bytes, SHA-256
`f32a3ffc325d7bca0511873a334047e964efe9a4b1e868ee50aba9fbedf1fdaa`.
The [independent reviewer](review_shared_threadripper_repaired_20261003.py)
passes runtime, environment, resources and native-output checks. Readback
verifies all four category pins and decisions. Frozen argv/input identities,
pre/post runtime bindings, raw resource counters and whole-run sampled
process/pressure observations agree with the retained evidence.

| Endpoint | Retained Observation |
| --- | ---: |
| Input proteins | 165,168 |
| Orthogroups and root HOGs, including singletons | 56,918 |
| Native ortholog-pair rows | 340,183 |
| Native wall seconds | 2,086.947282256 |
| Native CPU seconds | 58,240.723374 |
| Native-step lifetime peak bytes | 8,269,803,520 |
| Preflight foreign CPU-core equivalents | 53.9905288738 |
| Maximum sampled foreign CPU-core equivalents | 54.8075973213 |

Native-output validity is not a new accuracy evaluation. Wall time is the
independently replayed native monotonic command interval. CPU includes the
native-task-subtree wrapper bracket; memory is native-step lifetime peak
including launcher, not process RSS or whole-job teardown memory. Process and
pressure failure maps are empty. Contention distortion remains unknown and
potentially method dependent; no overhead subtraction or isolated ranking is
claimed. Preserve this repeat without fastest-repeat selection.

## Retained Repeat Difference

Repeat 0 ([22399](threadripper_shared_attempt_22399.json), index 3) has 56,919
groups, not 56,918. A bounded diagnostic compares unordered gene-membership
sets after dropping group labels. Exactly two repeat-0 groups are absent from
repeat 1: `{ENSMODP00000038126, ENSRNOP00000073598}` and
`{ENSRNOP00000060626}`. Exactly one repeat-1 group is absent from repeat 0:
`{ENSMODP00000038126, ENSRNOP00000060626, ENSRNOP00000073598}`.
All other membership sets agree. Both partitions contain all 165,168 proteins.

Parse the native pair tables with `csv.DictReader(delimiter='\t')`, normalize
each `(gene_a, gene_b)` to an unordered pair, and compare sets: both contain
340,183 unique pairs, with zero differences in either direction. Their raw
files are also byte-identical: 35,583,418 bytes, SHA-256
`05dcd749878009196dc9c6f8940650cb8e258e52fb35580af70c5c0414300328`.
The differing group files retain hashes
`87fed25ad8d3a283506e74bd11ec282a6119da358d672644acfcb16255083c04`
(repeat 0) and
`ed37210b635d6d9fb5f1116fab04d6d99cc26cb18112de061f08b58e703e5111`
(repeat 1). These identities are retained in the two native-output receipts.
This documents grouping variability, not its cause, deterministic grouping,
prediction accuracy or a reason to retune/retry the frozen timing panel.
Group co-membership must not be substituted for native ortholog pairs.

## Panel Continuation

The [fifteen-attempt table](threadripper_shared_panel_snapshot_20261003_v15/panel.json)
and [updated partial figure](threadripper_shared_resource_figure_20261003_v14/shared_threadripper_resources.pdf)
retain fourteen eligible observations and excluded index 0. All eight-proteome
cells now have two eligible repeats; no cell has three eligible repeats and
all medians/ranges remain unavailable. Earlier snapshots, reporting archive
and manuscript PDF retain their unchanged dated payloads.

Only after canonical terminal review does the launcher release index 15,
OrthoHMM high sensitivity on twelve proteomes, repeat 1, as job 22411. The
[live checkpoint](threadripper_shared_live_22411.json) is dated evidence, not a
terminal review; query that exact job on resumption. Review its terminal
evidence before index 16 (OrthoHMM satellite_v2, twelve proteomes, repeat 1).
Full-panel reporting, manuscript reconciliation and versioned release remain
unfinished. No scientific tuning, retry, quiet-window/DGX requirement or
unrelated workload change is introduced.
