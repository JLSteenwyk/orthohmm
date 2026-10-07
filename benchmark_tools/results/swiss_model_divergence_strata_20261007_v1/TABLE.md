# SwissTrees Fixed-Model Divergence Strata

Descriptive, development-exposed; WAG+G4 model-estimated distances, not biological time.
Scores are percentages; differences are percentage points. No subgroup confidence intervals.
Macro-family precision/recall then harmonic F1; initial HMM search on, profile refinement off.

Median family-distance cutoff: 2.2458246300500004

| Stratum | Families | Cell | F1 (%) | Precision (%) | Recall (%) |
| --- | ---: | --- | ---: | ---: | ---: |
| all | 18 | p0_c0_r0 | 68.918 | 64.394 | 74.127 |
| all | 18 | p0_c0_r1 | 78.957 | 94.915 | 67.593 |
| all | 18 | p0_c1_r0 | 68.571 | 60.871 | 78.502 |
| lower_or_equal_median | 9 | p0_c0_r0 | 69.563 | 61.605 | 79.883 |
| lower_or_equal_median | 9 | p0_c0_r1 | 81.508 | 94.022 | 71.934 |
| lower_or_equal_median | 9 | p0_c1_r0 | 68.318 | 57.058 | 85.116 |
| higher_than_median | 9 | p0_c0_r0 | 67.771 | 67.183 | 68.370 |
| higher_than_median | 9 | p0_c0_r1 | 76.199 | 95.808 | 63.253 |
| higher_than_median | 9 | p0_c1_r0 | 68.095 | 64.683 | 71.887 |

| Stratum | Families | Conditional Contrast | Delta F1 (pp) | Delta Precision (pp) | Delta Recall (pp) |
| --- | ---: | --- | ---: | ---: | ---: |
| all | 18 | R_at_P0_C0 | +10.039 | +30.521 | -6.533 |
| all | 18 | C_at_P0_R0 | -0.347 | -3.523 | +4.375 |
| lower_or_equal_median | 9 | R_at_P0_C0 | +11.945 | +32.417 | -7.950 |
| lower_or_equal_median | 9 | C_at_P0_R0 | -1.245 | -4.547 | +5.233 |
| higher_than_median | 9 | R_at_P0_C0 | +8.427 | +28.625 | -5.117 |
| higher_than_median | 9 | C_at_P0_R0 | +0.324 | -2.500 | +3.517 |

## Limitations

- Development-exposed descriptive strata, not independent confirmation or causal explanation.
- Fixed WAG+G4 and retained alignments; model fit and tree-search uncertainty remain unvalidated.
- Distance includes paralogy/sampling/composition/domain/alignment effects, not biological time.
- Macro-family precision/recall then harmonic F1, not pooled pairs or mean-family F1.
- Two conditional contrasts only; candidate-by-reconciliation interaction is not identified.
- No subgroup intervals, significance, new bootstrap draws, tuning or default change.
- Initial HMM search is on and downstream profile refinement is off in these three cells.
- Cell7's failed timing remains ineligible; new construction cost is not inference timing.
- Shared-host timing distortion is unknown and potentially tool-dependent.
