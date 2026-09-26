# Descriptive Duplication-Annotation Strata

Scores are percentages; differences are percentage points versus full OrthoFinder.
F1 is the harmonic mean of macro precision and recall, not mean family F1.
No new intervals or significance tests. Missing is not zero.

## all (18 families)

| Method | F1 | Precision | Recall | F1 difference | Status |
|---|---:|---:|---:|---:|---|
| orthohmm_high_sensitivity | 68.550 | 64.100 | 73.664 | -16.292 | descriptive |
| orthohmm_phylogeny_satellite_v2 | 83.351 | 95.518 | 73.934 | -1.490 | descriptive |
| orthofinder_3_1_5_full | 84.841 | 93.785 | 77.455 | 0.000 | descriptive |
| orthofinder_3_1_5_sequence_only | 69.051 | 56.907 | 87.786 | -15.790 | descriptive |
| sonicparanoid_2_0_9 | 79.846 | 87.228 | 73.616 | -4.995 | descriptive |
| proteinortho_6_3_6 | 71.811 | 95.218 | 57.641 | -13.030 | descriptive |
| fastoma_0_3_5 | 78.022 | 92.360 | 67.537 | -6.819 | descriptive |
| orthomcl_1_4 | 76.661 | 87.413 | 68.264 | -8.180 | descriptive |

## lower_duplication_fraction (9 families)

| Method | F1 | Precision | Recall | F1 difference | Status |
|---|---:|---:|---:|---:|---|
| orthohmm_high_sensitivity | 69.876 | 64.960 | 75.596 | -20.577 | descriptive |
| orthohmm_phylogeny_satellite_v2 | 86.378 | 97.943 | 77.256 | -4.075 | descriptive |
| orthofinder_3_1_5_full | 90.453 | 96.251 | 85.314 | 0.000 | descriptive |
| orthofinder_3_1_5_sequence_only | 69.349 | 56.399 | 90.020 | -21.103 | descriptive |
| sonicparanoid_2_0_9 | 82.314 | 88.273 | 77.108 | -8.139 | descriptive |
| proteinortho_6_3_6 | 76.145 | 97.410 | 62.501 | -14.308 | descriptive |
| fastoma_0_3_5 | 81.299 | 97.134 | 69.904 | -9.154 | descriptive |
| orthomcl_1_4 | 80.258 | 87.939 | 73.810 | -10.195 | descriptive |

## upper_duplication_fraction (9 families)

| Method | F1 | Precision | Recall | F1 difference | Status |
|---|---:|---:|---:|---:|---|
| orthohmm_high_sensitivity | 67.218 | 63.239 | 71.731 | -11.773 | descriptive |
| orthohmm_phylogeny_satellite_v2 | 80.309 | 93.092 | 70.613 | 1.318 | descriptive |
| orthofinder_3_1_5_full | 78.991 | 91.320 | 69.596 | 0.000 | descriptive |
| orthofinder_3_1_5_sequence_only | 68.715 | 57.415 | 85.552 | -10.277 | descriptive |
| sonicparanoid_2_0_9 | 77.328 | 86.183 | 70.124 | -1.663 | descriptive |
| proteinortho_6_3_6 | 67.350 | 93.026 | 52.782 | -11.641 | descriptive |
| fastoma_0_3_5 | 74.734 | 87.586 | 65.171 | -4.257 | descriptive |
| orthomcl_1_4 | 72.850 | 86.887 | 62.718 | -6.141 | descriptive |

## missing_duplication_fraction (0 families)

| Method | F1 | Precision | Recall | F1 difference | Status |
|---|---:|---:|---:|---:|---|
| orthohmm_high_sensitivity | NA | NA | NA | NA | empty_bin |
| orthohmm_phylogeny_satellite_v2 | NA | NA | NA | NA | empty_bin |
| orthofinder_3_1_5_full | NA | NA | NA | NA | empty_bin |
| orthofinder_3_1_5_sequence_only | NA | NA | NA | NA | empty_bin |
| sonicparanoid_2_0_9 | NA | NA | NA | NA | empty_bin |
| proteinortho_6_3_6 | NA | NA | NA | NA | empty_bin |
| fastoma_0_3_5 | NA | NA | NA | NA | empty_bin |
| orthomcl_1_4 | NA | NA | NA | NA | empty_bin |

Reference-derived, development-exposed annotation fractions are not evolutionary duplication rates or causal explanations. Default-S nodes are not explicit speciation observations. FastOMA uses a supplied tree; sequence-only OrthoFinder is a group-clique diagnostic.
