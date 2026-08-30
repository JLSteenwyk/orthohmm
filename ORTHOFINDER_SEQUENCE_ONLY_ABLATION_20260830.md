# OrthoFinder Sequence-Only Ablation, 2026-08-30

## Definition

OrthoFinder 3.1.5 does not expose a CLI option that exits after initial
orthogroup clustering. This ablation therefore uses the exact preserved
`clusters_OrthoFinder_I1.2.txt_id_pairs.txt` checkpoint from each completed
3.1.5 DIAMOND run. The checkpoint contains OrthoFinder's similarity graph
after MCL clustering but before multiple-sequence alignments, gene trees,
species-tree inference, STRIDE rooting, reconciliation, or hierarchical
orthogroup inference. No external species tree was supplied.

This is a controlled pre-phylogenetic ablation: sequence search, graph
construction, MCL parameters, and input proteins are identical to the full
OrthoFinder run. A validated converter maps the internal IDs through
`SequenceIDs.txt` and checks that every protein occurs in exactly one cluster.

## Primary Results

All metrics are higher-is-better. QfO is the mean of the six retained official
QfO 2020 assessments.

| method | OrthoBench F | precision | recall | exact | QfO mean |
| --- | ---: | ---: | ---: | ---: | ---: |
| OrthoFinder 3.1.5, sequence-only MCL | 58.7 | 45.7 | 82.0 | 18 | 0.544959 |
| OrthoHMM high sensitivity | 70.4 | 78.9 | 63.5 | 13 | 0.682548 |
| OrthoFinder 3.1.5, full phylogenetic | 72.7 | 66.1 | 80.9 | 19 | 0.782071 |

OrthoHMM beats sequence-only OrthoFinder by 11.7 OrthoBench F-score points
and 0.137589 QfO mean. Full OrthoFinder remains better than OrthoHMM by 2.3
OrthoBench F-score points and 0.099522 QfO mean.

## QfO Metrics

| metric | sequence-only OF | OrthoHMM | OrthoHMM delta | full OF |
| --- | ---: | ---: | ---: | ---: |
| VGNC F | 0.140300 | 0.667448 | +0.527148 | 0.988229 |
| SwissTrees F | 0.699507 | 0.673851 | -0.025657 | 0.859160 |
| TreeFam-A F | 0.711013 | 0.575870 | -0.135143 | 0.742887 |
| EC | 0.753394 | 0.931045 | +0.177650 | 0.941982 |
| GO | 0.407600 | 0.472282 | +0.064682 | 0.468166 |
| FAS | 0.557941 | 0.774794 | +0.216853 | 0.692001 |
| mean | 0.544959 | 0.682548 | +0.137589 | 0.782071 |

OrthoHMM wins four of six metrics and the mean. Sequence-only OrthoFinder
wins SwissTrees and TreeFam-A, but its very low VGNC precision and weaker
functional scores dominate the aggregate result.

## Interpretation

The pre-phylogenetic MCL groups maximize coverage and recall but over-group
aggressively. On OrthoBench, full OrthoFinder gains 20.4 precision points while
losing only 1.1 recall points relative to the MCL checkpoint, raising F-score
by 14.0 points. This directly quantifies the benefit of OrthoFinder's
phylogenetic stages on this dataset.

The QfO projection tells the same story. Sequence-only OrthoFinder emits
162,991,772 filtered cross-species pairs, 12.2 times full OrthoFinder's
13,403,651 and 18.7 times OrthoHMM's 8,720,518. The relation explosion also
made official scoring substantially more expensive: FAS used 45.6 GB peak RSS
and 58 minutes 40 seconds; the complete assessment took 1 hour 21 minutes.

For OrthoHMM development, permissive sequence clustering alone is not the
route to beating full OrthoFinder. OrthoHMM already has much better precision
than the sequence-only MCL result, but it remains too fragmented and has lower
recall than full OrthoFinder. The useful target is HMM-supported recovery of
remote homologs combined with duplication-aware safeguards that avoid the
MCL checkpoint's over-grouping.

## Reproducibility

The converter and its tests are in
`benchmark_tools/orthofinder_mcl_to_orthogroups.py` and
`tests/unit/test_orthofinder_mcl_to_orthogroups.py`. The official QfO scoring
job was 20868 and exited `0:0`; all six metrics completed. The unit suite
passes with 277 tests. Exact scores, paths, checksums, resource measurements,
and deltas are recorded in
`benchmark_tools/results/orthofinder_sequence_only_ablation_20260830.json`.

Inference runtime is not reported for the ablation because the MCL checkpoint
was extracted from the full run rather than produced by a separately timed
stop-after-MCL execution. Conversion and assessment times must not be treated
as sequence-only inference runtime.
