# Corrected Profile-Refined Reconciliation Scores

Scoring21783 completed0:0 in29:54 and independent admission21784
completed0:0 in13seconds. A fresh execution of the frozen admission script
reproduces the [retained receipt](qfo_corrected_factorial_score_admission_21784.json)
byte-for-byte, SHA-256
`a6c015fd7745c31f2c4fd25a70364bb126c2a55183306358e0d5bae34ca0c97b`.
The repeated receipt is retained at
`benchmarks/work/qfo_corrected_factorial_score_admission_recheck_21784.json`.
These elapsed times are shared-host workflow records, not comparative timings.

The p1_c0_r1 cell uses multi-sequence HMM profile refinement and inferred
reconciliation without candidate expansion. All5113180native phylogenetic
pairs map to the reference, with zero mapping losses and no root-HOG clique
expansion. It is not the full p1_c1_r1 publication configuration.

The fresh [seven-cell table](qfo_corrected_factorial_scores_20260919_v4/scores.md)
preserves the final missing cell and all earlier exports. Manifest SHA-256:
`5788b94ca021ee87b3d97d9bb68c2e989a411a2745bbad9ea0dc78a8209a1c13`.

## Profile Refinement With C-Off and R-On

| Endpoint | p0_c0_r1 | p1_c0_r1 |
| --- | ---: | ---: |
| GO similarity | 0.490260 | 0.490229 |
| EC similarity | 0.967702 | 0.968271 |
| VGNC F1 | 0.898185 | 0.897844 |
| SwissTrees F1 | 0.789574 | 0.786372 |
| TreeFam-A F1 | 0.602508 | 0.602432 |
| FAS | 0.785212 | 0.784576 |
| Project-defined secondary mean | 0.755573 | 0.754954 |

For p1_c0_r1, precision/recall are0.999385/0.815033 for VGNC,
0.949051/0.671303 for SwissTrees and0.956628/0.439650 for TreeFam-A.
Five of six point estimates decrease slightly and EC increases. These
development-exposed end-to-end contrasts do not establish statistical
significance, a direct mechanism, a default change or general superiority.
P-off retains initial HMM search; this is not an HMM-versus-no-HMM contrast.
GO/EC/FAS are not F1; the six-metric mean is a secondary project summary.

Complete-factorial uncertainty21894 still awaits final score admission21788.
The existing bootstrap seed, sample count, multiplicity adjustment and
endpoints remain unchanged; historical-input intervals are not transferred.
