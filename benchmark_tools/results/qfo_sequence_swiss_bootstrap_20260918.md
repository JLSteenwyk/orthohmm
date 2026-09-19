# Corrected QfO Sequence-Control SwissTrees Intervals

100000 shared family draws; seed 20260923; six adjusted endpoints.
Differences are candidate minus initial HMM, in raw 0-to-1 units.

| Contrast | Metric | Difference | Nominal 95% CI | Adjusted CI | Wins/ties/losses |
|---|---|---:|---|---|---|
| all_hits - p0_c0_r0 | F1 | -0.061092 | [-0.150733, 0.008469] | [-0.188023, 0.028017] | 5/4/9 |
| all_hits - p0_c0_r0 | PPV | -0.082110 | [-0.155468, -0.019045] | [-0.181914, -0.001387] | 7/4/7 |
| all_hits - p0_c0_r0 | TPR | -0.029194 | [-0.185829, 0.091725] | [-0.247769, 0.114154] | 11/4/3 |
| top100 - p0_c0_r0 | F1 | -0.059995 | [-0.149571, 0.009560] | [-0.186513, 0.029034] | 5/4/9 |
| top100 - p0_c0_r0 | PPV | -0.082123 | [-0.155473, -0.019056] | [-0.181927, -0.001404] | 7/4/7 |
| top100 - p0_c0_r0 | TPR | -0.026345 | [-0.182879, 0.093424] | [-0.245355, 0.114662] | 11/4/3 |

- Development-exposed exploratory analysis, not independent confirmation.
- Only 18 families; disjoint genes do not establish exchangeability.
- Adjustment covers these six endpoints only, not other QfO scores or their custom mean.
- Equal E-values do not establish matched sensitivity or cost.
- Source-bound counts were reconstructed; upstream inference/scoring admissions were not rerun.
