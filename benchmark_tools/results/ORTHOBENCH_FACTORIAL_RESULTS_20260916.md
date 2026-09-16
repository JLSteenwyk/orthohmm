# OrthoBench Factorial Uncertainty

Development-exposed component analysis, not independent confirmation.

70 RefOGs; 20,000 shared paired draws; seed 20260918.

| Cell | F1 (%) | Precision (%) | Recall (%) |
| --- | ---: | ---: | ---: |
| p0_c0_r0 | 69.763388 | 78.868592 | 62.542944 |
| p0_c0_r1 | 72.704996 | 87.494842 | 62.192223 |
| p0_c1_r0 | 66.639457 | 64.488265 | 68.939121 |
| p0_c1_r1 | 73.402288 | 81.486169 | 66.777581 |
| p1_c0_r0 | 70.358998 | 78.949544 | 63.454479 |
| p1_c0_r1 | 73.311364 | 87.458585 | 63.103758 |
| p1_c1_r0 | 67.207424 | 64.700133 | 69.916876 |
| p1_c1_r1 | 74.106074 | 81.770454 | 67.755336 |

Conditional effects are factor on minus factor off, in percentage points.

| Factor | On - Off | Metric | Difference | Nominal CI | Adjusted CI |
| --- | --- | --- | ---: | --- | --- |
| profile_expansion | p1_c0_r0 - p0_c0_r0 | f_score | 0.596 | [-0.750, 2.766] | [-1.164, 4.401] |
| profile_expansion | p1_c0_r0 - p0_c0_r0 | precision | 0.081 | [-0.739, 0.969] | [-1.193, 1.546] |
| profile_expansion | p1_c0_r0 - p0_c0_r0 | recall | 0.912 | [-0.776, 3.893] | [-1.184, 6.317] |
| profile_expansion | p1_c0_r1 - p0_c0_r1 | f_score | 0.606 | [-0.782, 2.829] | [-1.211, 4.521] |
| profile_expansion | p1_c0_r1 - p0_c0_r1 | precision | -0.036 | [-0.813, 0.668] | [-1.316, 1.145] |
| profile_expansion | p1_c0_r1 - p0_c0_r1 | recall | 0.912 | [-0.776, 3.893] | [-1.184, 6.317] |
| profile_expansion | p1_c1_r0 - p0_c1_r0 | f_score | 0.568 | [-0.695, 2.651] | [-1.113, 4.291] |
| profile_expansion | p1_c1_r0 - p0_c1_r0 | precision | 0.212 | [-0.697, 1.419] | [-1.151, 2.252] |
| profile_expansion | p1_c1_r0 - p0_c1_r0 | recall | 0.978 | [-0.745, 4.003] | [-1.158, 6.492] |
| profile_expansion | p1_c1_r1 - p0_c1_r1 | f_score | 0.704 | [-0.647, 2.907] | [-1.035, 4.579] |
| profile_expansion | p1_c1_r1 - p0_c1_r1 | precision | 0.284 | [-0.580, 1.289] | [-0.989, 1.986] |
| profile_expansion | p1_c1_r1 - p0_c1_r1 | recall | 0.978 | [-0.745, 4.003] | [-1.158, 6.492] |
| candidate_expansion | p0_c1_r0 - p0_c0_r0 | f_score | -3.124 | [-7.568, 0.896] | [-10.659, 3.275] |
| candidate_expansion | p0_c1_r0 - p0_c0_r0 | precision | -14.380 | [-21.259, -8.490] | [-26.493, -5.793] |
| candidate_expansion | p0_c1_r0 - p0_c0_r0 | recall | 6.396 | [2.885, 10.103] | [1.406, 12.573] |
| candidate_expansion | p0_c1_r1 - p0_c0_r1 | f_score | 0.697 | [-1.396, 3.037] | [-2.459, 4.821] |
| candidate_expansion | p0_c1_r1 - p0_c0_r1 | precision | -6.009 | [-9.387, -2.890] | [-11.610, -1.227] |
| candidate_expansion | p0_c1_r1 - p0_c0_r1 | recall | 4.585 | [2.186, 7.258] | [1.251, 9.273] |
| candidate_expansion | p1_c1_r0 - p1_c0_r0 | f_score | -3.152 | [-7.573, 0.899] | [-10.641, 3.249] |
| candidate_expansion | p1_c1_r0 - p1_c0_r0 | precision | -14.249 | [-21.156, -8.362] | [-26.433, -5.702] |
| candidate_expansion | p1_c1_r0 - p1_c0_r0 | recall | 6.462 | [2.931, 10.192] | [1.444, 12.708] |
| candidate_expansion | p1_c1_r1 - p1_c0_r1 | f_score | 0.795 | [-1.253, 3.108] | [-2.353, 4.857] |
| candidate_expansion | p1_c1_r1 - p1_c0_r1 | precision | -5.688 | [-9.019, -2.635] | [-11.290, -1.090] |
| candidate_expansion | p1_c1_r1 - p1_c0_r1 | recall | 4.652 | [2.228, 7.340] | [1.260, 9.358] |
| reconciliation | p0_c0_r1 - p0_c0_r0 | f_score | 2.942 | [0.841, 5.949] | [0.254, 8.144] |
| reconciliation | p0_c0_r1 - p0_c0_r0 | precision | 8.626 | [2.732, 16.019] | [0.786, 21.388] |
| reconciliation | p0_c0_r1 - p0_c0_r0 | recall | -0.351 | [-0.847, 0.000] | [-1.214, 0.000] |
| reconciliation | p0_c1_r1 - p0_c1_r0 | f_score | 6.763 | [3.159, 11.285] | [1.376, 14.612] |
| reconciliation | p0_c1_r1 - p0_c1_r0 | precision | 16.998 | [10.124, 24.781] | [6.890, 30.900] |
| reconciliation | p0_c1_r1 - p0_c1_r0 | recall | -2.162 | [-3.976, -0.498] | [-5.227, -0.062] |
| reconciliation | p1_c0_r1 - p1_c0_r0 | f_score | 2.952 | [0.845, 5.979] | [0.259, 8.158] |
| reconciliation | p1_c0_r1 - p1_c0_r0 | precision | 8.509 | [2.689, 15.853] | [0.769, 21.053] |
| reconciliation | p1_c0_r1 - p1_c0_r0 | recall | -0.351 | [-0.847, 0.000] | [-1.214, 0.000] |
| reconciliation | p1_c1_r1 - p1_c1_r0 | f_score | 6.899 | [3.271, 11.456] | [1.432, 14.808] |
| reconciliation | p1_c1_r1 - p1_c1_r0 | precision | 17.070 | [10.135, 24.881] | [6.872, 31.157] |
| reconciliation | p1_c1_r1 - p1_c1_r0 | recall | -2.162 | [-3.976, -0.498] | [-5.227, -0.062] |

Bonferroni adjustment retains all 36 prespecified contrast/metric endpoints.

## Descriptive Family Counts

| On - Off | Family F1 Wins | Ties | Losses |
| --- | ---: | ---: | ---: |
| p1_c0_r0 - p0_c0_r0 | 4 | 60 | 6 |
| p1_c0_r1 - p0_c0_r1 | 4 | 61 | 5 |
| p1_c1_r0 - p0_c1_r0 | 6 | 59 | 5 |
| p1_c1_r1 - p0_c1_r1 | 5 | 59 | 6 |
| p0_c1_r0 - p0_c0_r0 | 17 | 25 | 28 |
| p0_c1_r1 - p0_c0_r1 | 19 | 30 | 21 |
| p1_c1_r0 - p1_c0_r0 | 17 | 25 | 28 |
| p1_c1_r1 - p1_c0_r1 | 18 | 31 | 21 |
| p0_c0_r1 - p0_c0_r0 | 11 | 59 | 0 |
| p0_c1_r1 - p0_c1_r0 | 23 | 42 | 5 |
| p1_c0_r1 - p1_c0_r0 | 11 | 59 | 0 |
| p1_c1_r1 - p1_c1_r0 | 24 | 41 | 5 |

## Descriptive Interactions

Difference of conditional effects; no interaction confidence intervals or significance claims.

| Factors | Fixed Setting | F1 | Precision | Recall |
| --- | --- | ---: | ---: | ---: |
| profile_expansion x candidate_expansion | reconciliation=0 | -0.028 | 0.131 | 0.066 |
| profile_expansion x candidate_expansion | reconciliation=1 | 0.097 | 0.321 | 0.066 |
| profile_expansion x reconciliation | candidate_expansion=0 | 0.011 | -0.117 | 0.000 |
| profile_expansion x reconciliation | candidate_expansion=1 | 0.136 | 0.072 | -0.000 |
| candidate_expansion x reconciliation | profile_expansion=0 | 3.821 | 8.372 | -1.811 |
| candidate_expansion x reconciliation | profile_expansion=1 | 3.946 | 8.561 | -1.811 |

## Failed Cells

None.

## Limitations

- Development-exposed evidence; intervals are not selection-adjusted independent confirmation.
- Resampling units are RefOGs, not gene pairs; shared history and merged predictions can violate exchangeability.
- All 36 planned contrast/metric endpoints remain in the Bonferroni adjustment even when cells fail.
- Family wins and interaction differences are descriptive, without additional inferential claims.
- Weighted F1 is recomputed within every draw, not averaged across family F1 values.
- Percentile intervals are approximate; failed cells have no imputed score.
- Profile-off retains the HMM-based initial search; these contrasts do not measure the total HMM contribution.
- Reconciliation contrasts also change the output from candidate co-membership to final root HOGs.
- Execution, native-output conversion, reference coverage and provenance require separate validation.

## Coverage And Incremental Resources

All assigned genes includes singletons; multispecies coverage is not orthology accuracy.
Times are shared-node cached reconciliation costs, not end-to-end efficiency comparisons.
NA for candidate-only cells means no separate cell cost was measured, not zero cost.

| Cell | Groups | Singleton Groups | Genes In Multispecies Groups | Wall (s) | Mean CPU Cores | Peak Tree RSS (GiB) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| p0_c0_r0 | 63245 | 43417 | 192906 | NA | NA | NA |
| p0_c0_r1 | 64925 | 43882 | 192142 | 1582.81 | 23.12 | 1.467 |
| p0_c1_r0 | 54745 | 38468 | 200074 | NA | NA | NA |
| p0_c1_r1 | 60092 | 40581 | 197378 | 1952.12 | 26.74 | 1.573 |
| p1_c0_r0 | 62885 | 43358 | 193292 | NA | NA | NA |
| p1_c0_r1 | 64616 | 43845 | 192526 | 1514.54 | 25.51 | 1.476 |
| p1_c1_r0 | 54445 | 38355 | 200646 | NA | NA | NA |
| p1_c1_r1 | 59770 | 40521 | 197868 | 1963.79 | 28.67 | 1.577 |

Original reconciliation batch failures are retained in the machine-readable scheduler records.
Recovered postflight verification errors are not relabeled as native inference failures.
