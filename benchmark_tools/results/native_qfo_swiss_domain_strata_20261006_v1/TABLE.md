# Native SwissTrees Domain Bins

Descriptive percentages; differences R1 minus R0 in percentage points.
No new intervals, significance or causal domain effects. Initial HMM search on in both cells.

| Bin | Families | R0 F1 | R1 F1 | F1 change | Precision change | Recall change |
|---|---:|---:|---:|---:|---:|---:|
| all | 18 | 68.918 | 78.957 | +10.039 | +30.521 | -6.533 |
| median_pfam_types_below_two | 12 | 71.006 | 75.897 | +4.891 | +24.288 | -9.169 |
| median_pfam_types_at_least_two | 6 | 63.986 | 84.932 | +20.946 | +42.987 | -1.262 |
| repeated_type_fraction_below_quarter | 15 | 69.545 | 79.749 | +10.204 | +29.147 | -5.181 |
| repeated_type_fraction_at_least_quarter | 3 | 65.595 | 74.731 | +9.136 | +37.390 | -13.296 |

- Retrospective fixed-bin projections, not a new inferential test or method tuning.
- Macro precision/recall then harmonic F1, not pair-pooled or mean-family F1.
- Pfam counts are annotations, not complete validated domain architectures or causal mechanisms.
- Repeated-type high bin has three families; bins overlap; no new or historical intervals attached.
- FAS reference annotations do not independently validate FAS; no fragment or domain-loss truth.
- Direct source/raw/selected-annotation checks, not full transitive admission.
- R1 failed timing remains ineligible. Shared-host postprocessing time is not inference cost;
- contention effects are unknown and potentially tool-dependent, not isolated speed evidence.
