# Corrected Expanded Reconciliation Scores

Scoring21779completed0:0 in30:02; independent score admission21780completed
0:0 in13seconds. These are shared-host workflow records, not comparative
runtime measurements. The unchanged [admission receipt](qfo_corrected_factorial_score_admission_21780.json)
has SHA-256 `3f3683b893397ad056f85b5161498f4c8afc00f2b33e1ba6324e355ea1e835c1`.
A fresh execution of the same frozen admission script reproduced that entire
receipt byte-for-byte, retained at
`benchmarks/work/qfo_corrected_factorial_score_admission_recheck_21780.json`.

All5,977,100native phylogenetically inferred pairs were retained with zero
mapping losses. The native output is not expanded into root-HOG cliques.
The existing exporter now admits six of eight cells in a fresh
[machine-generated table](qfo_corrected_factorial_scores_20260919_v3/scores.md).
Its manifest SHA-256 is
`28d85568729cda1815d6bf06dafee08c80893a1ea17f04c97aa335a9c59b154a`.
Earlier partial tables remain unchanged; absent cells remain missing, not zero.

## Candidate Expansion With Reconciliation Enabled

Both cells retain initial HMM search and disable multi-sequence profile
expansion (P-off). Neither is the complete satellite_v2 publication method.

| Endpoint | p0_c0_r1 | p0_c1_r1 |
| --- | ---: | ---: |
| GO similarity | 0.490260 | 0.490035 |
| EC similarity | 0.967702 | 0.965435 |
| VGNC F1 | 0.898185 | 0.901626 |
| SwissTrees F1 | 0.789574 | 0.836354 |
| TreeFam-A F1 | 0.602508 | 0.616864 |
| FAS | 0.785212 | 0.762823 |
| Project-defined secondary mean | 0.755573 | 0.762189 |

| Reference | C-off precision | C-on precision | C-off recall | C-on recall |
| --- | ---: | ---: | ---: | ---: |
| VGNC | 0.999539 | 0.999491 | 0.815493 | 0.821217 |
| SwissTrees | 0.949152 | 0.955781 | 0.675932 | 0.743456 |
| TreeFam-A | 0.956682 | 0.959568 | 0.439719 | 0.454531 |

The observed SwissTrees F1 gain is4.677930percentage points. All three
reference-orthology F1 values increase, whereas GO,EC andFAS decrease. These
point estimates do not establish statistical significance, a full-factorial
interaction, a default change or general superiority. GO/EC/FAS are not F1,
and the six-endpoint mean is a project-defined secondary summary.

As already documented in the [native tree diagnostic](QFO_CORRECTED_EXPANDED_RECONCILIATION_20260919.md),
candidate expansion also changes the inferred species tree. This is an
end-to-end contrast, not isolation under a fixed realized tree; neither tree
is truth and topology change is not a demonstrated causal explanation.
The frozen complete-factorial uncertainty job21894awaits the final two
independent score admissions. Historical-release intervals are not transferred.
