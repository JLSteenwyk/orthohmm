# Stratified OrthoBench Error Analysis

Exploratory, development-exposed comparisons; descriptors do not establish biological mechanisms.

Both OrthoHMM configurations minus full OrthoFinder3.1.5. Differences in percentage points.
20,000 paired RefOG resamples; seed20260918; adjustment retains all84 planned endpoints.
Bins with fewer than five families have no bootstrap intervals. Empty bins are not scored as zero.

| Stratum | Families | Method | F1 (%) | F1 Difference | Nominal 95% CI | Adjusted CI | Status |
| --- | ---: | --- | ---: | ---: | --- | --- | --- |
| composition:concentrated | 1 | orthohmm_high_sensitivity | 22.388 | -35.762 | NA | NA | descriptive_only_lt5 |
| composition:concentrated | 1 | orthohmm_phylogeny_satellite_v2 | 36.622 | -21.528 | NA | NA | descriptive_only_lt5 |
| composition:missing | 0 | orthohmm_high_sensitivity | NA | NA | NA | NA | empty_nonestimable |
| composition:missing | 0 | orthohmm_phylogeny_satellite_v2 | NA | NA | NA | NA | empty_nonestimable |
| composition:not_concentrated | 69 | orthohmm_high_sensitivity | 72.716 | -1.270 | [-7.090, 4.453] | [-10.703, 8.774] | paired_bootstrap |
| composition:not_concentrated | 69 | orthohmm_phylogeny_satellite_v2 | 76.246 | 2.260 | [-3.861, 9.161] | [-7.418, 14.454] | paired_bootstrap |
| copy_number:multi_copy | 64 | orthohmm_high_sensitivity | 70.038 | -2.429 | [-8.356, 3.606] | [-12.775, 7.990] | paired_bootstrap |
| copy_number:multi_copy | 64 | orthohmm_phylogeny_satellite_v2 | 73.548 | 1.081 | [-4.965, 7.936] | [-8.452, 13.710] | paired_bootstrap |
| copy_number:single_copy | 6 | orthohmm_high_sensitivity | 78.978 | -1.649 | [-6.838, 4.929] | [-9.807, 7.548] | paired_bootstrap |
| copy_number:single_copy | 6 | orthohmm_phylogeny_satellite_v2 | 91.611 | 10.984 | [0.000, 17.908] | [0.000, 18.825] | paired_bootstrap |
| identity:higher_identity | 29 | orthohmm_high_sensitivity | 80.601 | -0.151 | [-9.586, 7.895] | [-16.484, 12.452] | paired_bootstrap |
| identity:higher_identity | 29 | orthohmm_phylogeny_satellite_v2 | 88.345 | 7.592 | [-2.969, 22.405] | [-6.407, 35.019] | paired_bootstrap |
| identity:lower_identity | 29 | orthohmm_high_sensitivity | 64.960 | -4.278 | [-14.583, 6.938] | [-20.913, 14.967] | paired_bootstrap |
| identity:lower_identity | 29 | orthohmm_phylogeny_satellite_v2 | 67.716 | -1.522 | [-10.657, 8.852] | [-16.103, 17.098] | paired_bootstrap |
| identity:missing | 12 | orthohmm_high_sensitivity | 65.106 | -4.601 | [-15.182, 5.759] | [-25.247, 14.735] | paired_bootstrap |
| identity:missing | 12 | orthohmm_phylogeny_satellite_v2 | 67.682 | -2.026 | [-10.760, 6.828] | [-19.553, 14.497] | paired_bootstrap |
| relative_length:missing | 0 | orthohmm_high_sensitivity | NA | NA | NA | NA | empty_nonestimable |
| relative_length:missing | 0 | orthohmm_phylogeny_satellite_v2 | NA | NA | NA | NA | empty_nonestimable |
| relative_length:not_short_relative | 30 | orthohmm_high_sensitivity | 72.446 | -2.337 | [-12.689, 8.272] | [-19.586, 15.687] | paired_bootstrap |
| relative_length:not_short_relative | 30 | orthohmm_phylogeny_satellite_v2 | 78.557 | 3.773 | [-4.052, 13.108] | [-8.231, 20.221] | paired_bootstrap |
| relative_length:short_relative | 40 | orthohmm_high_sensitivity | 69.174 | -2.468 | [-9.142, 4.123] | [-13.846, 9.317] | paired_bootstrap |
| relative_length:short_relative | 40 | orthohmm_phylogeny_satellite_v2 | 71.679 | 0.036 | [-7.783, 9.047] | [-12.824, 16.334] | paired_bootstrap |
| size:large_gt_50 | 6 | orthohmm_high_sensitivity | 51.483 | -16.469 | [-30.554, -5.004] | [-37.462, 0.291] | paired_bootstrap |
| size:large_gt_50 | 6 | orthohmm_phylogeny_satellite_v2 | 56.336 | -11.616 | [-23.430, -3.526] | [-33.873, 1.847] | paired_bootstrap |
| size:medium_21_50 | 30 | orthohmm_high_sensitivity | 77.289 | 1.232 | [-6.586, 9.135] | [-12.347, 15.501] | paired_bootstrap |
| size:medium_21_50 | 30 | orthohmm_phylogeny_satellite_v2 | 77.377 | 1.320 | [-5.948, 8.985] | [-10.497, 15.089] | paired_bootstrap |
| size:small_2_20 | 34 | orthohmm_high_sensitivity | 70.530 | -0.663 | [-12.423, 8.490] | [-20.916, 13.892] | paired_bootstrap |
| size:small_2_20 | 34 | orthohmm_phylogeny_satellite_v2 | 82.325 | 11.132 | [-3.133, 25.143] | [-8.550, 35.421] | paired_bootstrap |

## All Endpoint Effects

| Stratum | Method | Metric | Difference | Nominal 95% CI | Adjusted CI |
| --- | --- | --- | ---: | --- | --- |
| composition:concentrated | orthohmm_high_sensitivity | f_score | -35.762 | NA | NA |
| composition:concentrated | orthohmm_high_sensitivity | precision | 4.505 | NA | NA |
| composition:concentrated | orthohmm_high_sensitivity | recall | -63.045 | NA | NA |
| composition:concentrated | orthohmm_phylogeny_satellite_v2 | f_score | -21.528 | NA | NA |
| composition:concentrated | orthohmm_phylogeny_satellite_v2 | precision | 8.548 | NA | NA |
| composition:concentrated | orthohmm_phylogeny_satellite_v2 | recall | -49.961 | NA | NA |
| composition:missing | orthohmm_high_sensitivity | f_score | NA | NA | NA |
| composition:missing | orthohmm_high_sensitivity | precision | NA | NA | NA |
| composition:missing | orthohmm_high_sensitivity | recall | NA | NA | NA |
| composition:missing | orthohmm_phylogeny_satellite_v2 | f_score | NA | NA | NA |
| composition:missing | orthohmm_phylogeny_satellite_v2 | precision | NA | NA | NA |
| composition:missing | orthohmm_phylogeny_satellite_v2 | recall | NA | NA | NA |
| composition:not_concentrated | orthohmm_high_sensitivity | f_score | -1.270 | [-7.090, 4.453] | [-10.703, 8.774] |
| composition:not_concentrated | orthohmm_high_sensitivity | precision | 11.626 | [2.738, 19.466] | [-3.430, 25.250] |
| composition:not_concentrated | orthohmm_high_sensitivity | recall | -14.232 | [-22.586, -6.579] | [-29.434, -1.448] |
| composition:not_concentrated | orthohmm_phylogeny_satellite_v2 | f_score | 2.260 | [-3.861, 9.161] | [-7.418, 14.454] |
| composition:not_concentrated | orthohmm_phylogeny_satellite_v2 | precision | 14.888 | [4.443, 25.278] | [-1.188, 33.147] |
| composition:not_concentrated | orthohmm_phylogeny_satellite_v2 | recall | -10.552 | [-18.984, -3.094] | [-26.978, 1.668] |
| copy_number:multi_copy | orthohmm_high_sensitivity | f_score | -2.429 | [-8.356, 3.606] | [-12.775, 7.990] |
| copy_number:multi_copy | orthohmm_high_sensitivity | precision | 13.176 | [4.260, 20.711] | [-3.445, 25.924] |
| copy_number:multi_copy | orthohmm_high_sensitivity | recall | -17.591 | [-27.255, -8.179] | [-34.699, -2.778] |
| copy_number:multi_copy | orthohmm_phylogeny_satellite_v2 | f_score | 1.081 | [-4.965, 7.936] | [-8.452, 13.710] |
| copy_number:multi_copy | orthohmm_phylogeny_satellite_v2 | precision | 15.232 | [5.185, 24.912] | [-1.075, 32.353] |
| copy_number:multi_copy | orthohmm_phylogeny_satellite_v2 | recall | -13.187 | [-22.285, -4.451] | [-29.647, 0.560] |
| copy_number:single_copy | orthohmm_high_sensitivity | f_score | -1.649 | [-6.838, 4.929] | [-9.807, 7.548] |
| copy_number:single_copy | orthohmm_high_sensitivity | precision | 5.719 | [-0.906, 22.115] | [-3.832, 33.845] |
| copy_number:single_copy | orthohmm_high_sensitivity | recall | -13.065 | [-23.715, -3.232] | [-30.701, 0.000] |
| copy_number:single_copy | orthohmm_phylogeny_satellite_v2 | f_score | 10.984 | [0.000, 17.908] | [0.000, 18.825] |
| copy_number:single_copy | orthohmm_phylogeny_satellite_v2 | precision | 30.788 | [0.000, 52.231] | [0.000, 61.491] |
| copy_number:single_copy | orthohmm_phylogeny_satellite_v2 | recall | -12.031 | [-28.244, 0.000] | [-41.935, 0.000] |
| identity:higher_identity | orthohmm_high_sensitivity | f_score | -0.151 | [-9.586, 7.895] | [-16.484, 12.452] |
| identity:higher_identity | orthohmm_high_sensitivity | precision | 3.574 | [-12.329, 13.372] | [-23.832, 17.326] |
| identity:higher_identity | orthohmm_high_sensitivity | recall | -5.409 | [-12.038, 1.578] | [-16.052, 7.605] |
| identity:higher_identity | orthohmm_phylogeny_satellite_v2 | f_score | 7.592 | [-2.969, 22.405] | [-6.407, 35.019] |
| identity:higher_identity | orthohmm_phylogeny_satellite_v2 | precision | 18.272 | [0.344, 39.720] | [-4.580, 55.219] |
| identity:higher_identity | orthohmm_phylogeny_satellite_v2 | recall | -5.047 | [-11.716, 1.830] | [-15.688, 7.365] |
| identity:lower_identity | orthohmm_high_sensitivity | f_score | -4.278 | [-14.583, 6.938] | [-20.913, 14.967] |
| identity:lower_identity | orthohmm_high_sensitivity | precision | 22.513 | [10.670, 30.523] | [1.166, 36.944] |
| identity:lower_identity | orthohmm_high_sensitivity | recall | -24.928 | [-39.553, -7.384] | [-47.694, 2.265] |
| identity:lower_identity | orthohmm_phylogeny_satellite_v2 | f_score | -1.522 | [-10.657, 8.852] | [-16.103, 17.098] |
| identity:lower_identity | orthohmm_phylogeny_satellite_v2 | precision | 17.919 | [4.432, 28.592] | [-4.792, 36.990] |
| identity:lower_identity | orthohmm_phylogeny_satellite_v2 | recall | -19.026 | [-32.915, -2.805] | [-41.933, 6.067] |
| identity:missing | orthohmm_high_sensitivity | f_score | -4.601 | [-15.182, 5.759] | [-25.247, 14.735] |
| identity:missing | orthohmm_high_sensitivity | precision | 8.456 | [-7.680, 30.394] | [-19.058, 42.306] |
| identity:missing | orthohmm_high_sensitivity | recall | -16.674 | [-30.751, -6.803] | [-48.932, -0.098] |
| identity:missing | orthohmm_phylogeny_satellite_v2 | f_score | -2.026 | [-10.760, 6.828] | [-19.553, 14.497] |
| identity:missing | orthohmm_phylogeny_satellite_v2 | precision | 6.286 | [-6.240, 24.376] | [-13.929, 35.454] |
| identity:missing | orthohmm_phylogeny_satellite_v2 | recall | -10.715 | [-26.157, -0.837] | [-43.433, 5.569] |
| relative_length:missing | orthohmm_high_sensitivity | f_score | NA | NA | NA |
| relative_length:missing | orthohmm_high_sensitivity | precision | NA | NA | NA |
| relative_length:missing | orthohmm_high_sensitivity | recall | NA | NA | NA |
| relative_length:missing | orthohmm_phylogeny_satellite_v2 | f_score | NA | NA | NA |
| relative_length:missing | orthohmm_phylogeny_satellite_v2 | precision | NA | NA | NA |
| relative_length:missing | orthohmm_phylogeny_satellite_v2 | recall | NA | NA | NA |
| relative_length:not_short_relative | orthohmm_high_sensitivity | f_score | -2.337 | [-12.689, 8.272] | [-19.586, 15.687] |
| relative_length:not_short_relative | orthohmm_high_sensitivity | precision | 10.629 | [-10.091, 22.220] | [-21.511, 27.522] |
| relative_length:not_short_relative | orthohmm_high_sensitivity | recall | -14.852 | [-32.875, 5.251] | [-42.165, 13.691] |
| relative_length:not_short_relative | orthohmm_phylogeny_satellite_v2 | f_score | 3.773 | [-4.052, 13.108] | [-8.231, 20.221] |
| relative_length:not_short_relative | orthohmm_phylogeny_satellite_v2 | precision | 18.678 | [0.900, 28.551] | [-4.329, 36.446] |
| relative_length:not_short_relative | orthohmm_phylogeny_satellite_v2 | recall | -10.132 | [-24.961, 5.775] | [-32.926, 13.585] |
| relative_length:short_relative | orthohmm_high_sensitivity | f_score | -2.468 | [-9.142, 4.123] | [-13.846, 9.317] |
| relative_length:short_relative | orthohmm_high_sensitivity | precision | 13.975 | [4.387, 23.191] | [-1.979, 29.796] |
| relative_length:short_relative | orthohmm_high_sensitivity | recall | -18.901 | [-29.099, -10.100] | [-36.475, -5.208] |
| relative_length:short_relative | orthohmm_phylogeny_satellite_v2 | f_score | 0.036 | [-7.783, 9.047] | [-12.824, 16.334] |
| relative_length:short_relative | orthohmm_phylogeny_satellite_v2 | precision | 14.074 | [1.674, 26.993] | [-3.840, 37.986] |
| relative_length:short_relative | orthohmm_phylogeny_satellite_v2 | recall | -14.833 | [-25.791, -5.456] | [-34.240, -0.370] |
| size:large_gt_50 | orthohmm_high_sensitivity | f_score | -16.469 | [-30.554, -5.004] | [-37.462, 0.291] |
| size:large_gt_50 | orthohmm_high_sensitivity | precision | 24.372 | [9.554, 33.936] | [2.751, 40.741] |
| size:large_gt_50 | orthohmm_high_sensitivity | recall | -39.931 | [-62.261, -15.623] | [-78.262, -11.260] |
| size:large_gt_50 | orthohmm_phylogeny_satellite_v2 | f_score | -11.616 | [-23.430, -3.526] | [-33.873, 1.847] |
| size:large_gt_50 | orthohmm_phylogeny_satellite_v2 | precision | 14.183 | [-1.527, 23.793] | [-4.876, 30.239] |
| size:large_gt_50 | orthohmm_phylogeny_satellite_v2 | recall | -31.728 | [-57.020, -6.064] | [-74.865, 1.898] |
| size:medium_21_50 | orthohmm_high_sensitivity | f_score | 1.232 | [-6.586, 9.135] | [-12.347, 15.501] |
| size:medium_21_50 | orthohmm_high_sensitivity | precision | 12.283 | [1.587, 24.072] | [-2.797, 31.629] |
| size:medium_21_50 | orthohmm_high_sensitivity | recall | -8.601 | [-18.113, 0.596] | [-24.977, 6.654] |
| size:medium_21_50 | orthohmm_phylogeny_satellite_v2 | f_score | 1.320 | [-5.948, 8.985] | [-10.497, 15.089] |
| size:medium_21_50 | orthohmm_phylogeny_satellite_v2 | precision | 8.611 | [-2.515, 21.066] | [-7.975, 30.240] |
| size:medium_21_50 | orthohmm_phylogeny_satellite_v2 | recall | -5.735 | [-14.202, 2.496] | [-20.440, 8.289] |
| size:small_2_20 | orthohmm_high_sensitivity | f_score | -0.663 | [-12.423, 8.490] | [-20.916, 13.892] |
| size:small_2_20 | orthohmm_high_sensitivity | precision | 6.879 | [-11.958, 18.626] | [-26.223, 25.995] |
| size:small_2_20 | orthohmm_high_sensitivity | recall | -13.547 | [-21.462, -5.712] | [-26.671, -0.608] |
| size:small_2_20 | orthohmm_phylogeny_satellite_v2 | f_score | 11.132 | [-3.133, 25.143] | [-8.550, 35.421] |
| size:small_2_20 | orthohmm_phylogeny_satellite_v2 | precision | 26.773 | [3.845, 44.128] | [-2.350, 56.494] |
| size:small_2_20 | orthohmm_phylogeny_satellite_v2 | recall | -10.144 | [-18.591, -1.964] | [-24.788, 3.581] |

## Limitations

- Post-development exploratory analysis; no independent confirmation, equivalence claim or new default.
- Strata overlap and share histories; within-stratum RefOG resampling assumes exchangeable families, not independent gene pairs.
- All84 endpoints remain in multiplicity adjustment, including empty and descriptive-only bins.
- Adjusted percentile tails have about six draws per tail with20,000 replicates; uncertainty is approximate.
- Copy number, identity, relative length and global composition do not establish duplication history, fragments, domains or causal mechanisms.
- Sufficient statistics retain full-reference scoring and low-certainty conventions before restriction to each stratum.
- This analysis uses previously audited retained group outputs, not new end-to-end inference or QfO evidence.
