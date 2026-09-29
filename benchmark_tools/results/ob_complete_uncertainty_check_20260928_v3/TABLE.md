# Complete OrthoBench Paired Comparisons

Exploratory, development-exposed; conditional on exchangeable RefOGs. Differences in percentage points versus full OrthoFinder 3.1.5.

| Method | Metric | Difference | Nominal 95% CI | Adjusted CI (21 endpoints) |
|---|---|---:|---|---|
| fastoma_0_3_5 | f_score | -41.830 | [-50.523, -32.079] | [-54.989, -25.918] |
| fastoma_0_3_5 | precision | 27.566 | [16.284, 36.862] | [10.202, 41.760] |
| fastoma_0_3_5 | recall | -62.398 | [-70.002, -53.667] | [-73.635, -48.383] |
| orthofinder_3_1_5_sequence_only | f_score | -14.031 | [-21.998, -5.493] | [-26.657, -1.859] |
| orthofinder_3_1_5_sequence_only | precision | -20.359 | [-30.871, -8.773] | [-37.537, -3.154] |
| orthofinder_3_1_5_sequence_only | recall | 1.133 | [0.000, 3.566] | [0.000, 5.546] |
| orthohmm_high_sensitivity | f_score | -2.377 | [-8.102, 3.453] | [-11.045, 6.746] |
| orthohmm_high_sensitivity | precision | 12.884 | [4.336, 20.275] | [-0.895, 24.222] |
| orthohmm_high_sensitivity | recall | -17.452 | [-26.889, -8.224] | [-32.103, -3.845] |
| orthohmm_phylogeny_satellite_v2 | f_score | 1.370 | [-4.482, 8.051] | [-7.281, 12.166] |
| orthohmm_phylogeny_satellite_v2 | precision | 15.705 | [5.994, 25.186] | [1.175, 30.634] |
| orthohmm_phylogeny_satellite_v2 | recall | -13.151 | [-21.995, -4.625] | [-26.947, -0.609] |
| orthomcl_1_4 | f_score | -17.671 | [-27.573, -6.980] | [-32.754, -1.040] |
| orthomcl_1_4 | precision | -6.987 | [-24.706, 15.076] | [-34.019, 25.711] |
| orthomcl_1_4 | recall | -29.344 | [-38.936, -19.468] | [-43.938, -14.781] |
| proteinortho_6_3_6 | f_score | -27.679 | [-35.610, -18.555] | [-39.688, -12.857] |
| proteinortho_6_3_6 | precision | 30.991 | [19.851, 40.162] | [13.704, 44.843] |
| proteinortho_6_3_6 | recall | -51.568 | [-59.065, -43.054] | [-62.477, -38.295] |
| sonicparanoid_2_0_9 | f_score | -25.979 | [-45.101, -0.028] | [-52.286, 5.451] |
| sonicparanoid_2_0_9 | precision | -29.530 | [-51.551, 15.826] | [-60.591, 25.358] |
| sonicparanoid_2_0_9 | recall | -15.984 | [-26.681, -5.935] | [-32.653, -2.246] |

100,000 paired draws; PCG64 seed 20260928. Approximate percentile intervals, not guarantees of coverage. No tuning or new independent confirmation.

Family F1 comparisons use rational counts and the retained 1e-10 percentage-point tie tolerance.

| Method | Family F1 wins | Ties | Losses |
|---|---:|---:|---:|
| fastoma_0_3_5 | 6 | 5 | 59 |
| orthofinder_3_1_5_sequence_only | 0 | 60 | 10 |
| orthohmm_high_sensitivity | 20 | 11 | 39 |
| orthohmm_phylogeny_satellite_v2 | 23 | 10 | 37 |
| orthomcl_1_4 | 15 | 13 | 42 |
| proteinortho_6_3_6 | 7 | 8 | 55 |
| sonicparanoid_2_0_9 | 19 | 19 | 32 |
