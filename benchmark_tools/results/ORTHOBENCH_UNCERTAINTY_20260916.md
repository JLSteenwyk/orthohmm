# OrthoBench Paired Uncertainty

Development-exposed analysis; not independent confirmation or a superiority claim.

Baseline: `orthofinder_3_1_5_full`. 70 RefOGs; 20,000 paired bootstrap replicates; seed 20260916.

| Method | F1 (%) | Precision (%) | Recall (%) |
| --- | ---: | ---: | ---: |
| orthofinder_3_1_5_full | 72.736480 | 66.065103 | 80.906577 |
| orthohmm_high_sensitivity | 70.358998 | 78.949544 | 63.454479 |
| orthohmm_phylogeny_satellite_v2 | 74.106074 | 81.770454 | 67.755336 |

All differences below are method minus baseline, in percentage points.

| Method | Metric | Difference | Paired 95% CI | Multiplicity-adjusted CI |
| --- | --- | ---: | --- | --- |
| orthohmm_high_sensitivity | f_score | -2.377 | [-8.178, 3.348] | [-10.135, 5.478] |
| orthohmm_high_sensitivity | precision | 12.884 | [4.211, 20.220] | [1.281, 22.627] |
| orthohmm_high_sensitivity | recall | -17.452 | [-26.965, -8.214] | [-30.347, -5.540] |
| orthohmm_phylogeny_satellite_v2 | f_score | 1.370 | [-4.504, 7.917] | [-6.290, 10.627] |
| orthohmm_phylogeny_satellite_v2 | precision | 15.705 | [5.983, 25.105] | [2.802, 28.663] |
| orthohmm_phylogeny_satellite_v2 | recall | -13.151 | [-22.064, -4.566] | [-25.248, -2.242] |

Bonferroni tail adjustment over 6 reported contrasts/metrics.

| Method | Family F1 wins | Ties | Losses |
| --- | ---: | ---: | ---: |
| orthohmm_high_sensitivity | 20 | 11 | 39 |
| orthohmm_phylogeny_satellite_v2 | 23 | 10 | 37 |

## Limitations

- Development-exposed benchmark; intervals do not correct for previous method selection.
- RefOG resampling assumes exchangeable families; shared histories and fused predictions can violate independence.
- Percentile intervals are approximate, not a guarantee of simultaneous coverage.
- Family F1 wins are descriptive; the benchmark is not a mean of family F1 values.
- No gene-pair independence assumption and no bootstrap-derived p-values are used.
