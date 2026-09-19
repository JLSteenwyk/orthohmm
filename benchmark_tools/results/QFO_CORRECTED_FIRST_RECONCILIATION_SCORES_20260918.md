# First Corrected Reconciliation Scores

Scoring job 21775 completed 0:0 in 30:12; independent admission 21776
completed 0:0 in 12 seconds. The unchanged
[admission receipt](qfo_corrected_factorial_score_admission_21776.json) has
SHA-256 `f36c1518cbdaa50160f837e8909efc8fb771f207fa66b684ad322c61fb5cfdcc`.
All 5,113,820 native phylogenetic pairs map without loss. These elapsed times
are shared-host workflow records, not comparative efficiency measurements.

The existing hash-checking exporter generated the
[updated eight-cell table](qfo_corrected_factorial_scores_20260918_v2/scores.md)
and manifest from the five admitted receipts. Three cells remain missing,
not zero. Version 1 remains unchanged.

| Endpoint | p0_c0_r0 | p0_c0_r1 |
| --- | ---: | ---: |
| GO similarity | 0.472119 | 0.490260 |
| EC similarity | 0.932114 | 0.967702 |
| VGNC F1 | 0.666834 | 0.898185 |
| SwissTrees F1 | 0.689184 | 0.789574 |
| TreeFam-A F1 | 0.605404 | 0.602508 |
| FAS | 0.777125 | 0.785212 |
| Project-defined secondary mean | 0.690463 | 0.755573 |

In this fixed P-off/C-off contrast, reconciliation increases precision and
reduces recall in each reference-orthology benchmark:

| Reference | R-off precision | R-on precision | R-off recall | R-on recall |
| --- | ---: | ---: | ---: | ---: |
| VGNC | 0.555120 | 0.999539 | 0.834837 | 0.815493 |
| SwissTrees | 0.643940 | 0.949152 | 0.741266 | 0.675932 |
| TreeFam-A | 0.792245 | 0.956682 | 0.489873 | 0.439719 |

These are admitted point estimates, not evidence of a statistically established
benefit. The prespecified full-factorial SwissTrees uncertainty analysis still
requires all eight cells. Original-release intervals cannot be transferred.
P-off retains initial HMM search; R-off predictions are group-derived clique
pairs whereas R-on predictions are native phylogenetic pairs. Neither arm
here is the full publication satellite_v2 configuration. No new parameter
choice or superiority claim follows from this development-exposed contrast.
