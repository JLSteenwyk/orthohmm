# Species-Tree Robustness

Exploratory prespecified topology stress test; unchanged frozen method and reused raw gene trees.

# OrthoBench Paired Uncertainty

Development-exposed analysis; not independent confirmation or a superiority claim.

Baseline: `supplied_control`. 70 RefOGs; 20,000 paired bootstrap replicates; seed 20260918.

| Method | F1 (%) | Precision (%) | Recall (%) |
| --- | ---: | ---: | ---: |
| nni1_0 | 74.236657 | 82.291583 | 67.618020 |
| nni1_1 | 74.031050 | 81.791309 | 67.615750 |
| nni1_2 | 74.029855 | 81.600959 | 67.744391 |
| nni2_0 | 74.015490 | 82.249855 | 67.279835 |
| nni2_1 | 74.270946 | 82.577826 | 67.482567 |
| nni2_2 | 73.744102 | 81.992367 | 67.003669 |
| supplied_control | 74.106074 | 81.770454 | 67.755336 |

All differences below are method minus baseline, in percentage points.

| Method | Metric | Difference | Paired 95% CI | Multiplicity-adjusted CI |
| --- | --- | ---: | --- | --- |
| nni1_0 | f_score | 0.131 | [-0.143, 0.538] | [-0.234, 0.823] |
| nni1_0 | precision | 0.521 | [-0.004, 1.535] | [-0.018, 2.361] |
| nni1_0 | recall | -0.137 | [-0.338, 0.000] | [-0.458, 0.000] |
| nni1_1 | f_score | -0.075 | [-0.250, 0.000] | [-0.390, 0.000] |
| nni1_1 | precision | 0.021 | [0.000, 0.089] | [-0.007, 0.158] |
| nni1_1 | recall | -0.140 | [-0.435, 0.000] | [-0.629, 0.000] |
| nni1_2 | f_score | -0.076 | [-1.382, 1.445] | [-2.077, 2.485] |
| nni1_2 | precision | -0.169 | [-3.812, 3.572] | [-5.733, 6.061] |
| nni1_2 | recall | -0.011 | [-0.523, 0.510] | [-0.836, 0.832] |
| nni2_0 | f_score | -0.091 | [-0.765, 0.478] | [-1.175, 0.772] |
| nni2_0 | precision | 0.479 | [-0.081, 1.507] | [-0.213, 2.317] |
| nni2_0 | recall | -0.476 | [-1.261, 0.000] | [-1.776, 0.000] |
| nni2_1 | f_score | 0.165 | [-0.200, 0.622] | [-0.342, 0.921] |
| nni2_1 | precision | 0.807 | [0.057, 1.938] | [-0.037, 2.810] |
| nni2_1 | recall | -0.273 | [-0.569, 0.000] | [-0.723, 0.000] |
| nni2_2 | f_score | -0.362 | [-0.958, 0.020] | [-1.347, 0.069] |
| nni2_2 | precision | 0.222 | [-0.020, 0.613] | [-0.069, 0.885] |
| nni2_2 | recall | -0.752 | [-1.762, 0.000] | [-2.375, 0.000] |

Bonferroni tail adjustment over 18 reported contrasts/metrics.

| Method | Family F1 wins | Ties | Losses |
| --- | ---: | ---: | ---: |
| nni1_0 | 2 | 66 | 2 |
| nni1_1 | 0 | 69 | 1 |
| nni1_2 | 3 | 62 | 5 |
| nni2_0 | 3 | 65 | 2 |
| nni2_1 | 4 | 63 | 3 |
| nni2_2 | 1 | 66 | 3 |

## Limitations

- Development-exposed benchmark; intervals do not correct for previous method selection.
- RefOG resampling assumes exchangeable families; shared histories and fused predictions can violate independence.
- Percentile intervals are approximate, not a guarantee of simultaneous coverage.
- Family F1 wins are descriptive; the benchmark is not a mean of family F1 values.
- No gene-pair independence assumption and no bootstrap-derived p-values are used.
- All six fixed variants and 18 F1/precision/recall endpoints are retained; no best-tree selection.
- The unchanged supplied-tree control exactly reproduces the inferred baseline; this is not independent tree reconstruction.
- Rooted NNI variants are controlled perturbations, not posterior draws or an empirical tree-error distribution.
- Lengths travel with subtrees; these topology edits do not preserve all evolutionary distances.
- Specified after development and YGOB outcomes; not independent confirmation, parameter robustness, or a new default.

## Coverage And Incremental Resources

All inputs include singletons. Costs reuse raw gene trees on a shared node; they are not end-to-end timings.

| Tree | Rooted Clade Distance | Groups | Singleton Groups | Genes In Multispecies Groups | Wall (s) | Peak Tree RSS (GiB) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| supplied_control | 0 | 59770 | 40521 | 197868 | 107.42 | 1.609 |
| nni1_0 | 2 | 59867 | 40554 | 197813 | 97.16 | 1.600 |
| nni1_1 | 2 | 59827 | 40541 | 197828 | 108.29 | 1.603 |
| nni1_2 | 2 | 58940 | 40204 | 198323 | 105.94 | 1.606 |
| nni2_0 | 4 | 59908 | 40579 | 197780 | 113.19 | 1.596 |
| nni2_1 | 4 | 60081 | 40686 | 197624 | 134.81 | 1.576 |
| nni2_2 | 4 | 59938 | 40591 | 197755 | 122.70 | 1.587 |
