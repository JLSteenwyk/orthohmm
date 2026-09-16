# Simulation Results: fixed_length_v1

Scores are percentages and means of completed seed-level metrics, not pooled gene pairs.
Failed or inapplicable seeds have no imputed score. Available-case means are conditional on success.

| Condition | Method | Completed | Failed | Inapplicable | Precision | Recall | F1 |
|---|---|---:|---:|---:|---:|---:|---:|
| baseline | OrthoHMM high sensitivity | 10 | 0 | 0 | 92.95 | 99.97 | 96.29 |
| baseline | OrthoHMM satellite_v2 | 10 | 0 | 0 | 99.99 | 99.37 | 99.68 |
| baseline | OrthoFinder full | 0 | 10 | 0 | NA | NA | NA |
| baseline | OrthoFinder sequence checkpoint | 0 | 10 | 0 | NA | NA | NA |
| divergent | OrthoHMM high sensitivity | 10 | 0 | 0 | 97.07 | 68.93 | 75.36 |
| divergent | OrthoHMM satellite_v2 | 7 | 3 | 0 | 99.85 | 89.88 | 94.09 |
| divergent | OrthoFinder full | 0 | 10 | 0 | NA | NA | NA |
| divergent | OrthoFinder sequence checkpoint | 0 | 10 | 0 | NA | NA | NA |
| turnover | OrthoHMM high sensitivity | 10 | 0 | 0 | 77.79 | 99.98 | 87.12 |
| turnover | OrthoHMM satellite_v2 | 10 | 0 | 0 | 98.65 | 96.73 | 97.66 |
| turnover | OrthoFinder full | 0 | 10 | 0 | NA | NA | NA |
| turnover | OrthoFinder sequence checkpoint | 0 | 10 | 0 | NA | NA | NA |
| divergent_turnover | OrthoHMM high sensitivity | 10 | 0 | 0 | 89.42 | 69.06 | 72.35 |
| divergent_turnover | OrthoHMM satellite_v2 | 7 | 3 | 0 | 99.53 | 89.37 | 93.67 |
| divergent_turnover | OrthoFinder full | 0 | 10 | 0 | NA | NA | NA |
| divergent_turnover | OrthoFinder sequence checkpoint | 0 | 10 | 0 | NA | NA | NA |
| missing20 | OrthoHMM high sensitivity | 10 | 0 | 0 | 93.48 | 99.94 | 96.56 |
| missing20 | OrthoHMM satellite_v2 | 10 | 0 | 0 | 99.81 | 99.26 | 99.53 |
| missing20 | OrthoFinder full | 0 | 10 | 0 | NA | NA | NA |
| missing20 | OrthoFinder sequence checkpoint | 0 | 10 | 0 | NA | NA | NA |
| uneven_taxa | OrthoHMM high sensitivity | 10 | 0 | 0 | 93.75 | 99.83 | 96.65 |
| uneven_taxa | OrthoHMM satellite_v2 | 10 | 0 | 0 | 99.91 | 98.71 | 99.30 |
| uneven_taxa | OrthoFinder full | 0 | 10 | 0 | NA | NA | NA |
| uneven_taxa | OrthoFinder sequence checkpoint | 0 | 10 | 0 | NA | NA | NA |
| taxon_count_control | OrthoHMM high sensitivity | 10 | 0 | 0 | 92.97 | 99.84 | 96.24 |
| taxon_count_control | OrthoHMM satellite_v2 | 10 | 0 | 0 | 99.97 | 98.35 | 99.15 |
| taxon_count_control | OrthoFinder full | 0 | 10 | 0 | NA | NA | NA |
| taxon_count_control | OrthoFinder sequence checkpoint | 0 | 10 | 0 | NA | NA | NA |

## Paired F1 Differences

OrthoHMM minus full OrthoFinder, in percentage points, using only successful paired seeds.

| Condition | Method | Paired Seeds | Difference | Nominal 95% CI | Bonferroni-14 CI |
|---|---|---:|---:|---|---|
| baseline | OrthoHMM high sensitivity | 0 | NA | NA | NA |
| baseline | OrthoHMM satellite_v2 | 0 | NA | NA | NA |
| divergent | OrthoHMM high sensitivity | 0 | NA | NA | NA |
| divergent | OrthoHMM satellite_v2 | 0 | NA | NA | NA |
| turnover | OrthoHMM high sensitivity | 0 | NA | NA | NA |
| turnover | OrthoHMM satellite_v2 | 0 | NA | NA | NA |
| divergent_turnover | OrthoHMM high sensitivity | 0 | NA | NA | NA |
| divergent_turnover | OrthoHMM satellite_v2 | 0 | NA | NA | NA |
| missing20 | OrthoHMM high sensitivity | 0 | NA | NA | NA |
| missing20 | OrthoHMM satellite_v2 | 0 | NA | NA | NA |
| uneven_taxa | OrthoHMM high sensitivity | 0 | NA | NA | NA |
| uneven_taxa | OrthoHMM satellite_v2 | 0 | NA | NA | NA |
| taxon_count_control | OrthoHMM high sensitivity | 0 | NA | NA | NA |
| taxon_count_control | OrthoHMM satellite_v2 | 0 | NA | NA | NA |

The sequence checkpoint requires a valid full parent run and has no independent runtime.
Intervals use 20,000 whole-seed bootstrap replicates; ten seeds limit tail resolution.
Resource logs are descriptive shared-machine measurements, not controlled scaling results.
This panel is not pooled with the other simulation panel and does not establish publication readiness.

Source results: `simulation_fixed_native_results_20260916.json`, SHA-256 `305dcb1dde0c0f57d8b390b0f00cc6148d95103a7efb98bd9f6dce37b96e64a8`.
