# OrthoBench Parameter Neighborhood

Post-development six-variant parameter sensitivity; no default selection.

# OrthoBench Paired Uncertainty

Development-exposed analysis; not independent confirmation or a superiority claim.

Baseline: `control`. 70 RefOGs; 20,000 paired bootstrap replicates; seed 20260918.

| Method | F1 (%) | Precision (%) | Recall (%) |
| --- | ---: | ---: | ---: |
| control | 74.106074 | 81.770454 | 67.755336 |
| cpm_high | 71.463468 | 78.271132 | 65.745248 |
| cpm_low | 74.586147 | 80.929105 | 69.165205 |
| margin_high | 73.274559 | 81.801415 | 66.357549 |
| margin_low | 74.973144 | 83.990775 | 67.704121 |
| norm_high | 74.110826 | 81.782028 | 67.755336 |
| norm_low | 74.008570 | 81.533400 | 67.755336 |

All differences below are method minus baseline, in percentage points.

| Method | Metric | Difference | Paired 95% CI | Multiplicity-adjusted CI |
| --- | --- | ---: | --- | --- |
| cpm_high | f_score | -2.643 | [-6.143, -0.167] | [-8.450, 0.708] |
| cpm_high | precision | -3.499 | [-10.865, 1.486] | [-16.021, 3.067] |
| cpm_high | recall | -2.010 | [-4.462, -0.437] | [-6.291, -0.148] |
| cpm_low | f_score | 0.480 | [-0.471, 1.548] | [-0.923, 2.180] |
| cpm_low | precision | -0.841 | [-2.025, 0.502] | [-2.665, 1.284] |
| cpm_low | recall | 1.410 | [0.218, 2.839] | [-0.139, 3.686] |
| margin_high | f_score | -0.832 | [-3.015, 0.451] | [-4.607, 0.792] |
| margin_high | precision | 0.031 | [-1.207, 1.392] | [-1.801, 2.329] |
| margin_high | recall | -1.398 | [-4.386, 0.000] | [-6.747, 0.000] |
| margin_low | f_score | 0.867 | [-0.702, 2.882] | [-1.362, 4.316] |
| margin_low | precision | 2.220 | [-1.964, 7.663] | [-3.392, 11.043] |
| margin_low | recall | -0.051 | [-0.721, 0.443] | [-1.175, 0.666] |
| norm_high | f_score | 0.005 | [0.000, 0.016] | [0.000, 0.024] |
| norm_high | precision | 0.012 | [0.000, 0.040] | [0.000, 0.062] |
| norm_high | recall | 0.000 | [0.000, 0.000] | [0.000, 0.000] |
| norm_low | f_score | -0.098 | [-0.230, 0.000] | [-0.318, 0.000] |
| norm_low | precision | -0.237 | [-0.552, 0.000] | [-0.771, 0.000] |
| norm_low | recall | 0.000 | [0.000, 0.000] | [0.000, 0.000] |

Bonferroni tail adjustment over 18 planned endpoints.

| Method | Family F1 wins | Ties | Losses |
| --- | ---: | ---: | ---: |
| cpm_high | 9 | 47 | 14 |
| cpm_low | 13 | 47 | 10 |
| margin_high | 9 | 54 | 7 |
| margin_low | 8 | 55 | 7 |
| norm_high | 1 | 69 | 0 |
| norm_low | 0 | 67 | 3 |

## Limitations

- Development-exposed benchmark; intervals do not correct for previous method selection.
- RefOG resampling assumes exchangeable families; shared histories and fused predictions can violate independence.
- Percentile intervals are approximate, not a guarantee of simultaneous coverage.
- Family F1 wins are descriptive; the benchmark is not a mean of family F1 values.
- No gene-pair independence assumption and no bootstrap-derived p-values are used.
- Post-development prespecified parameter sensitivity, not independent confirmation or default selection.
- All six variants and 18 endpoints remain in the planned family even when an inference run fails.
- Failed runs have no accuracy estimate or confidence interval; they are not assigned zero scores.
- Intervals including zero do not demonstrate equivalence or general robustness.
- CPM variants recompute their own HMM seeds and candidates; species trees are inferred under unchanged native rules. Reported phylogeny costs omit upstream replay and candidate construction and are not matched timing comparisons.


Incremental inferred-phylogeny runs with exact-input checkpoint reuse on a shared node; not end-to-end.
