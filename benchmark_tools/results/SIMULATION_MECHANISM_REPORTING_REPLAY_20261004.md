# Simulation Mechanism Summary

Finite-panel reporting replay from checked counts; not raw admission, native inference, confidence intervals or independent validation.

| Condition | Inferred F1 (%) | Generating-root F1 (%) | Change (pp) | Oracle FN | Cross-candidate FN | Different graph components | Connected but separated |
|---|---:|---:|---:|---:|---:|---:|---:|
| baseline | 99.493 | 99.699 | +0.206 | 177 | 177 | 177 | 0 |
| divergent | 69.631 | 69.768 | +0.136 | 12,766 | 12,736 | 11,442 | 1,294 |
| divergent_turnover | 71.672 | 72.168 | +0.496 | 13,104 | 13,072 | 11,554 | 1,518 |
| missing20 | 99.403 | 99.622 | +0.219 | 138 | 138 | 138 | 0 |
| taxon_count_control | 99.316 | 99.554 | +0.238 | 92 | 92 | 92 | 0 |
| turnover | 98.990 | 99.724 | +0.734 | 196 | 196 | 196 | 0 |
| uneven_taxa | 99.028 | 99.203 | +0.174 | 162 | 162 | 162 | 0 |

## Complete Within-Candidate Residual Cohort

| Mechanism | Error Pairs |
|---|---:|
| single_copy_bypass_on_true_duplication | 45 |
| true_duplication_without_retained_species_overlap | 21 |
| unsupported_satellite_constraint | 62 |

Screened 70 cells / 10,125 candidates; 1,377 eligible tree controls. Upstream: 163,527 true-pair rows.
Residual cohort: 20 candidates / 800 pair rows.
Missing significant hits do not isolate search-stage causes; candidate paths may cross other families.
Retained candidate membership and constraint policy are fixed. No defaults or original scores change.
