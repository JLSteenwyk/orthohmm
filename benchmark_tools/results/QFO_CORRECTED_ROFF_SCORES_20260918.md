# Corrected QfO Without Reconciliation

Four independently audited corrected-release assessments are admitted.
The retained reports are byte-for-byte copies of the admission outputs;
the export checks their hashes, conversion bindings and endpoint arithmetic.

| Cell | Scoring job | Admission job | SwissTrees F1 | Secondary mean | Mapped pairs |
| --- | --- | --- | ---: | ---: | ---: |
| p0_c0_r0 | 21773 | 21774 | 0.689184 | 0.690463 | 9,009,082 |
| p0_c1_r0 | 21777 | 21778 | 0.685709 | 0.684563 | 11,734,021 |
| p1_c0_r0 | 21781 | 21782 | 0.685498 | 0.689555 | 9,032,719 |
| p1_c1_r0 | 21785 | 21786 | 0.686087 | 0.684765 | 11,755,521 |

All four admission jobs completed with exit 0:0. All predictions are
cross-species group-derived clique pairs, with zero mapping losses.
P denotes multi-sequence profile expansion, C candidate-family expansion,
and R phylogenetic reconciliation. P-off still uses initial HMM search.

[Machine-generated six-endpoint table](qfo_corrected_factorial_scores_20260918_v1/scores.md)
and its manifest retain complete precision/recall details, source hashes and
prediction semantics. GO/EC similarity and FAS are not F1. The six-metric mean
is project-defined and secondary. Pair volume is not protein coverage.

Reconciled cells remain missing, not zero. This partial table cannot establish
candidate-expansion/reconciliation interactions. Historical-input scores are
not substituted, and no parameter selection is authorized by these results.
Shared-host elapsed times are not matched comparative efficiency evidence.

The p0_c0_r0 admission enables the already frozen three-arm SwissTrees
sequence-control analysis: initial HMM, DIAMOND all-hit, and DIAMOND top100,
with identical downstream settings. Full factorial uncertainty still requires
all eight admitted cells.
