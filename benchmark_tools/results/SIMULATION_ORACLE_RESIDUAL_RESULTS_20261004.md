# Generating-Root Residual Error Trace

All within-candidate generating-root FP/FN candidates in all 70 retained cells; post hoc diagnostic

Screened 70 cells / 10,125 candidates.

| Cell | Candidate | Status | Genes | TP | FP | FN | Error Classes |
|---|---|---|---:|---:|---:|---:|---|
| divergent_20261102 | Family0000133 | oracle_eligible | 16 | 41 | 0 | 15 | unsupported_satellite_constraint: 15 |
| divergent_20261105 | Family0000200 | oracle_eligible | 12 | 19 | 0 | 15 | unsupported_satellite_constraint: 15 |
| divergent_turnover_20261102 | Family0000108 | unambiguous_bypass | 8 | 13 | 15 | 0 | single_copy_bypass_on_true_duplication: 15 |
| divergent_turnover_20261102 | Family0000148 | oracle_eligible | 16 | 41 | 0 | 15 | unsupported_satellite_constraint: 15 |
| divergent_turnover_20261102 | Family0000192 | unambiguous_bypass | 6 | 10 | 5 | 0 | single_copy_bypass_on_true_duplication: 5 |
| divergent_turnover_20261105 | Family0000013 | oracle_eligible | 14 | 44 | 0 | 2 | unsupported_satellite_constraint: 2 |
| divergent_turnover_20261106 | Family0000128 | oracle_eligible | 9 | 21 | 14 | 0 | true_duplication_without_retained_species_overlap: 14 |
| divergent_turnover_20261107 | Family0000080 | unambiguous_bypass | 5 | 6 | 4 | 0 | single_copy_bypass_on_true_duplication: 4 |
| divergent_turnover_20261108 | Family0000013 | oracle_eligible | 5 | 5 | 4 | 0 | true_duplication_without_retained_species_overlap: 4 |
| divergent_turnover_20261108 | Family0000014 | oracle_eligible | 4 | 2 | 3 | 0 | true_duplication_without_retained_species_overlap: 3 |
| divergent_turnover_20261108 | Family0000126 | unambiguous_bypass | 4 | 4 | 2 | 0 | single_copy_bypass_on_true_duplication: 2 |
| divergent_turnover_20261108 | Family0000416 | unambiguous_bypass | 3 | 1 | 2 | 0 | single_copy_bypass_on_true_duplication: 2 |
| divergent_turnover_20261110 | Family0000067 | unambiguous_bypass | 2 | 0 | 1 | 0 | single_copy_bypass_on_true_duplication: 1 |
| divergent_turnover_20261110 | Family0000138 | oracle_eligible | 23 | 91 | 0 | 15 | unsupported_satellite_constraint: 15 |
| divergent_turnover_20261110 | Family0000142 | unambiguous_bypass | 5 | 4 | 6 | 0 | single_copy_bypass_on_true_duplication: 6 |
| divergent_turnover_20261110 | Family0000286 | unambiguous_bypass | 6 | 10 | 5 | 0 | single_copy_bypass_on_true_duplication: 5 |
| divergent_turnover_20261110 | Family0000306 | unambiguous_bypass | 3 | 1 | 2 | 0 | single_copy_bypass_on_true_duplication: 2 |
| missing20_20261103 | Family0000053 | unambiguous_bypass | 8 | 27 | 1 | 0 | single_copy_bypass_on_true_duplication: 1 |
| missing20_20261103 | Family0000085 | unambiguous_bypass | 7 | 20 | 1 | 0 | single_copy_bypass_on_true_duplication: 1 |
| missing20_20261108 | Family0000080 | unambiguous_bypass | 7 | 20 | 1 | 0 | single_copy_bypass_on_true_duplication: 1 |

## Aggregate Error Counts

| Mechanism Class | Pairs |
|---|---:|
| single_copy_bypass_on_true_duplication | 45 |
| true_duplication_without_retained_species_overlap | 21 |
| unsupported_satellite_constraint | 62 |

All cross-species pairs in selected candidates are traced, including correctly classified pairs.
Counts are descriptive simulator-specific mechanism evidence, not independent accuracy gains or a tuning prescription.
