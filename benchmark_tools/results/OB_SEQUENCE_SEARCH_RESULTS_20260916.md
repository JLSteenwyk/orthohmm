# Sequence-Search Control

Exploratory initial-search replacement; profiles, candidate expansion and reconciliation off.

# OrthoBench Paired Uncertainty

Development-exposed analysis; not independent confirmation or a superiority claim.

Baseline: `p0_c0_r0`. 70 RefOGs; 20,000 paired bootstrap replicates; seed 20260918.

| Method | F1 (%) | Precision (%) | Recall (%) |
| --- | ---: | ---: | ---: |
| all_hits | 65.762183 | 55.036231 | 81.680884 |
| p0_c0_r0 | 69.763388 | 78.868592 | 62.542944 |
| top100 | 66.484803 | 55.897861 | 82.019046 |

All differences below are method minus baseline, in percentage points.

| Method | Metric | Difference | Paired 95% CI | Multiplicity-adjusted CI |
| --- | --- | ---: | --- | --- |
| all_hits | f_score | -4.001 | [-12.050, 4.005] | [-14.834, 6.688] |
| all_hits | precision | -23.832 | [-32.033, -14.844] | [-34.605, -11.751] |
| all_hits | recall | 19.138 | [11.962, 26.386] | [9.289, 28.920] |
| top100 | f_score | -3.279 | [-11.528, 5.015] | [-14.424, 7.844] |
| top100 | precision | -22.971 | [-31.275, -13.723] | [-33.959, -10.442] |
| top100 | recall | 19.476 | [12.298, 26.694] | [9.600, 29.229] |

Bonferroni tail adjustment over 6 reported contrasts/metrics.

| Method | Family F1 wins | Ties | Losses |
| --- | ---: | ---: | ---: |
| all_hits | 39 | 10 | 21 |
| top100 | 39 | 10 | 21 |

## Limitations

- Development-exposed benchmark; intervals do not correct for previous method selection.
- RefOG resampling assumes exchangeable families; shared histories and fused predictions can violate independence.
- Percentile intervals are approximate, not a guarantee of simultaneous coverage.
- Family F1 wins are descriptive; the benchmark is not a mean of family F1 values.
- No gene-pair independence assumption and no bootstrap-derived p-values are used.
- All-hit DIAMOND is the primary control; post-search top100 is a reporting diagnostic, not HMM prefilter emulation.
- Equal E-value cutoffs and length divisors do not establish matched sensitivity or score calibration.
- Costs are shared-node incremental graph replays plus separately recorded search; not matched efficiency evidence.
- Control specified after development and YGOB outcomes; not independent confirmation or a new default.
