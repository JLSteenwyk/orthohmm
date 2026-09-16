# Simulation Results: variable_length_v2

Scores are percentages and means of completed seed-level metrics, not pooled gene pairs.
Failed or inapplicable seeds have no imputed score. Available-case means are conditional on success.

| Condition | Method | Completed | Failed | Inapplicable | Precision | Recall | F1 |
|---|---|---:|---:|---:|---:|---:|---:|
| baseline | OrthoHMM high sensitivity | 10 | 0 | 0 | 94.92 | 99.37 | 97.08 |
| baseline | OrthoHMM satellite_v2 | 10 | 0 | 0 | 99.92 | 99.02 | 99.46 |
| baseline | OrthoFinder full | 10 | 0 | 0 | 99.94 | 99.95 | 99.94 |
| baseline | OrthoFinder sequence checkpoint | 10 | 0 | 0 | 95.71 | 100.00 | 97.79 |
| divergent | OrthoHMM high sensitivity | 10 | 0 | 0 | 97.93 | 52.81 | 65.91 |
| divergent | OrthoHMM satellite_v2 | 9 | 1 | 0 | 99.84 | 57.30 | 70.64 |
| divergent | OrthoFinder full | 5 | 5 | 0 | 99.74 | 81.79 | 87.86 |
| divergent | OrthoFinder sequence checkpoint | 5 | 5 | 0 | 97.18 | 82.02 | 86.91 |
| turnover | OrthoHMM high sensitivity | 10 | 0 | 0 | 79.96 | 99.46 | 88.30 |
| turnover | OrthoHMM satellite_v2 | 10 | 0 | 0 | 99.13 | 98.56 | 98.84 |
| turnover | OrthoFinder full | 10 | 0 | 0 | 99.33 | 99.40 | 99.36 |
| turnover | OrthoFinder sequence checkpoint | 10 | 0 | 0 | 85.43 | 99.60 | 91.89 |
| divergent_turnover | OrthoHMM high sensitivity | 10 | 0 | 0 | 90.76 | 52.98 | 64.59 |
| divergent_turnover | OrthoHMM satellite_v2 | 8 | 2 | 0 | 98.87 | 61.55 | 74.34 |
| divergent_turnover | OrthoFinder full | 10 | 0 | 0 | 98.74 | 76.55 | 84.45 |
| divergent_turnover | OrthoFinder sequence checkpoint | 10 | 0 | 0 | 89.81 | 77.21 | 81.21 |
| missing20 | OrthoHMM high sensitivity | 10 | 0 | 0 | 94.85 | 99.25 | 96.98 |
| missing20 | OrthoHMM satellite_v2 | 10 | 0 | 0 | 99.85 | 98.96 | 99.40 |
| missing20 | OrthoFinder full | 10 | 0 | 0 | 99.45 | 99.68 | 99.56 |
| missing20 | OrthoFinder sequence checkpoint | 10 | 0 | 0 | 95.53 | 99.96 | 97.68 |
| uneven_taxa | OrthoHMM high sensitivity | 10 | 0 | 0 | 94.69 | 98.51 | 96.54 |
| uneven_taxa | OrthoHMM satellite_v2 | 10 | 0 | 0 | 100.00 | 98.07 | 98.98 |
| uneven_taxa | OrthoFinder full | 10 | 0 | 0 | 99.87 | 99.68 | 99.78 |
| uneven_taxa | OrthoFinder sequence checkpoint | 10 | 0 | 0 | 96.14 | 99.92 | 97.98 |
| taxon_count_control | OrthoHMM high sensitivity | 10 | 0 | 0 | 95.04 | 99.14 | 97.02 |
| taxon_count_control | OrthoHMM satellite_v2 | 10 | 0 | 0 | 99.90 | 98.35 | 99.10 |
| taxon_count_control | OrthoFinder full | 10 | 0 | 0 | 99.88 | 99.79 | 99.83 |
| taxon_count_control | OrthoFinder sequence checkpoint | 10 | 0 | 0 | 95.87 | 100.00 | 97.87 |

## Paired F1 Differences

OrthoHMM minus full OrthoFinder, in percentage points, using only successful paired seeds.

| Condition | Method | Paired Seeds | Difference | Nominal 95% CI | Bonferroni-14 CI |
|---|---|---:|---:|---|---|
| baseline | OrthoHMM high sensitivity | 10 | -2.87 | [-4.45, -1.54] | [-5.38, -1.14] |
| baseline | OrthoHMM satellite_v2 | 10 | -0.48 | [-1.12, -0.05] | [-1.51, -0.01] |
| divergent | OrthoHMM high sensitivity | 5 | -14.29 | [-20.75, -9.08] | [-22.46, -6.87] |
| divergent | OrthoHMM satellite_v2 | 5 | -11.74 | [-18.44, -5.04] | [-21.07, -3.24] |
| turnover | OrthoHMM high sensitivity | 10 | -11.06 | [-15.42, -7.17] | [-17.76, -5.74] |
| turnover | OrthoHMM satellite_v2 | 10 | -0.52 | [-0.96, -0.14] | [-1.18, 0.01] |
| divergent_turnover | OrthoHMM high sensitivity | 10 | -19.86 | [-23.77, -15.54] | [-25.38, -13.58] |
| divergent_turnover | OrthoHMM satellite_v2 | 8 | -12.28 | [-15.85, -8.47] | [-17.23, -6.68] |
| missing20 | OrthoHMM high sensitivity | 10 | -2.58 | [-4.11, -1.33] | [-5.03, -0.93] |
| missing20 | OrthoHMM satellite_v2 | 10 | -0.17 | [-0.83, 0.36] | [-1.21, 0.55] |
| uneven_taxa | OrthoHMM high sensitivity | 10 | -3.24 | [-5.73, -1.53] | [-7.25, -1.14] |
| uneven_taxa | OrthoHMM satellite_v2 | 10 | -0.80 | [-2.37, 0.12] | [-3.44, 0.21] |
| taxon_count_control | OrthoHMM high sensitivity | 10 | -2.81 | [-4.54, -1.43] | [-5.62, -1.04] |
| taxon_count_control | OrthoHMM satellite_v2 | 10 | -0.73 | [-1.75, -0.06] | [-2.48, -0.01] |

The sequence checkpoint requires a valid full parent run and has no independent runtime.
Intervals use 20,000 whole-seed bootstrap replicates; ten seeds limit tail resolution.
Resource logs are descriptive shared-machine measurements, not controlled scaling results.
This panel is not pooled with the other simulation panel and does not establish publication readiness.

Source results: `simulation_variable_native_results_20260916.json`, SHA-256 `cc99fc31c3433809098212d3c8dc12f4ad829b4cfeb66dabcb13aede4e95258f`.
