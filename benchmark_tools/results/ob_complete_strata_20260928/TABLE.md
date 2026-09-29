# OrthoBench All-Method Strata

Weighted F1 percentages; development-exposed descriptive results. No new confidence intervals.

Method columns: M1 = orthofinder_3_1_5_full; M2 = orthohmm_high_sensitivity; M3 = orthohmm_phylogeny_satellite_v2; M4 = orthofinder_3_1_5_sequence_only; M5 = sonicparanoid_2_0_9; M6 = proteinortho_6_3_6; M7 = fastoma_0_3_5; M8 = orthomcl_1_4

| Stratum | Families | M1 | M2 | M3 | M4 | M5 | M6 | M7 | M8 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| composition:concentrated | 1 | 58.150 | 22.388 | 36.622 | 58.150 | 8.190 | 5.061 | 5.063 | 6.965 |
| composition:missing | 0 | NA | NA | NA | NA | NA | NA | NA | NA |
| composition:not_concentrated | 69 | 73.986 | 72.716 | 76.246 | 58.743 | 47.756 | 47.282 | 32.469 | 56.917 |
| copy_number:multi_copy | 64 | 72.468 | 70.038 | 73.548 | 58.116 | 46.015 | 42.867 | 27.676 | 53.698 |
| copy_number:single_copy | 6 | 80.627 | 78.978 | 91.611 | 80.627 | 70.677 | 92.813 | 92.885 | 94.737 |
| identity:higher_identity | 29 | 80.752 | 80.601 | 88.345 | 61.713 | 73.196 | 66.941 | 48.203 | 63.545 |
| identity:lower_identity | 29 | 69.238 | 64.960 | 67.716 | 53.965 | 30.824 | 34.907 | 23.779 | 56.447 |
| identity:missing | 12 | 69.708 | 65.106 | 67.682 | 67.646 | 70.652 | 31.632 | 20.013 | 37.876 |
| relative_length:missing | 0 | NA | NA | NA | NA | NA | NA | NA | NA |
| relative_length:not_short_relative | 30 | 74.783 | 72.446 | 78.557 | 62.699 | 30.104 | 54.497 | 39.189 | 67.564 |
| relative_length:short_relative | 40 | 71.642 | 69.174 | 71.679 | 56.626 | 70.234 | 39.238 | 25.905 | 49.758 |
| size:large_gt_50 | 6 | 67.952 | 51.483 | 56.336 | 51.117 | 46.059 | 15.055 | 9.393 | 33.106 |
| size:medium_21_50 | 30 | 76.057 | 77.289 | 77.377 | 62.062 | 38.188 | 36.867 | 15.796 | 53.386 |
| size:small_2_20 | 34 | 71.193 | 70.530 | 82.325 | 60.190 | 80.296 | 78.467 | 68.791 | 72.475 |

Precision, recall, full-precision F1 and empty-bin status are in scores.tsv.
Family descriptors are proxies, not validated domain, fragment or duplication-history labels.
