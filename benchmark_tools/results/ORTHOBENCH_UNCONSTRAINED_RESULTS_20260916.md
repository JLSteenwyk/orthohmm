# Unconstrained Reconciliation Diagnostic

Exploratory unconstrained diagnostic; not a ninth factorial cell.

# OrthoBench Paired Uncertainty

Development-exposed analysis; not independent confirmation or a superiority claim.

Baseline: `p1_c1_r1`. 70 RefOGs; 20,000 paired bootstrap replicates; seed 20260918.

| Method | F1 (%) | Precision (%) | Recall (%) |
| --- | ---: | ---: | ---: |
| p1_c1_r1 | 74.106074 | 81.770454 | 67.755336 |
| p1_c1_r1_unconstrained_v2 | 72.614510 | 76.092405 | 69.440641 |

All differences below are method minus baseline, in percentage points.

| Method | Metric | Difference | Paired 95% CI | Multiplicity-adjusted CI |
| --- | --- | ---: | --- | --- |
| p1_c1_r1_unconstrained_v2 | f_score | -1.492 | [-4.841, 0.644] | [-5.873, 0.931] |
| p1_c1_r1_unconstrained_v2 | precision | -5.678 | [-12.424, -1.217] | [-14.307, -0.794] |
| p1_c1_r1_unconstrained_v2 | recall | 1.685 | [0.060, 3.409] | [0.013, 3.804] |

Bonferroni tail adjustment over 3 reported contrasts/metrics.

| Method | Family F1 wins | Ties | Losses |
| --- | ---: | ---: | ---: |
| p1_c1_r1_unconstrained_v2 | 5 | 55 | 10 |

## Limitations

- Development-exposed benchmark; intervals do not correct for previous method selection.
- RefOG resampling assumes exchangeable families; shared histories and fused predictions can violate independence.
- Percentile intervals are approximate, not a guarantee of simultaneous coverage.
- Family F1 wins are descriptive; the benchmark is not a mean of family F1 values.
- No gene-pair independence assumption and no bootstrap-derived p-values are used.
- Detailed execution specified after factorial outcomes; not independent confirmation or a new default.
