# Corrected SwissTrees Identity Strata

Descriptive percentages; differences are percentage points versus full OrthoFinder.
F1 is the harmonic mean of macro precision and recall. No additional intervals or significance claims.

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

## higher_identity (9 families)

| Method | F1 | Precision | Recall | F1 difference | Status |
|---|---:|---:|---:|---:|---|
| orthohmm_high_sensitivity | 72.435 | 63.932 | 83.547 | -20.165 | descriptive |
| orthohmm_phylogeny_satellite_v2 | 91.914 | 98.195 | 86.389 | -0.685 | descriptive |
| orthofinder_3_1_5_full | 92.600 | 97.697 | 88.007 | 0.000 | descriptive |
| orthofinder_3_1_5_sequence_only | 72.561 | 58.696 | 95.002 | -20.038 | descriptive |
| sonicparanoid_2_0_9 | 86.055 | 90.288 | 82.201 | -6.544 | descriptive |
| proteinortho_6_3_6 | 81.853 | 97.796 | 70.379 | -10.747 | descriptive |
| fastoma_0_3_5 | 85.377 | 92.810 | 79.046 | -7.223 | descriptive |
| orthomcl_1_4 | 83.426 | 83.513 | 83.340 | -9.173 | descriptive |

## lower_identity (9 families)

| Method | F1 | Precision | Recall | F1 difference | Status |
|---|---:|---:|---:|---:|---|
| orthohmm_high_sensitivity | 64.023 | 64.268 | 63.780 | -12.682 | descriptive |
| orthohmm_phylogeny_satellite_v2 | 73.973 | 92.840 | 61.480 | -2.732 | descriptive |
| orthofinder_3_1_5_full | 76.705 | 89.873 | 66.902 | 0.000 | descriptive |
| orthofinder_3_1_5_sequence_only | 65.457 | 55.117 | 80.571 | -11.248 | descriptive |
| sonicparanoid_2_0_9 | 73.372 | 84.167 | 65.031 | -3.333 | descriptive |
| proteinortho_6_3_6 | 60.488 | 92.640 | 44.904 | -16.217 | descriptive |
| fastoma_0_3_5 | 69.618 | 91.910 | 56.029 | -7.087 | descriptive |
| orthomcl_1_4 | 67.222 | 91.313 | 53.189 | -9.483 | descriptive |

## missing_identity (0 families)

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

Development-exposed, alignment-dependent identity bins are not calibrated evolutionary distances. Missing is not zero. FastOMA uses a supplied tree; sequence-only OrthoFinder is a group-clique diagnostic.
