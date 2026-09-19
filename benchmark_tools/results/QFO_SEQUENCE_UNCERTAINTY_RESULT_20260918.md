# Corrected QfO Sequence-Control Uncertainty

The frozen production analysis completed successfully for all 18 SwissTrees
families and 10,765 reference relations. Counts were reconstructed from
source-bound raw evidence for corrected initial HMM p0_c0_r0, DIAMOND
all-hit and DIAMOND top100 with the frozen downstream grouping settings.
No families were dropped and no historical prediction counts substituted.

| Search | F1 | Precision | Recall |
| --- | ---: | ---: | ---: |
| Initial HMM | 0.689184 | 0.643940 | 0.741266 |
| DIAMOND all-hit | 0.628092 | 0.561830 | 0.712073 |
| DIAMOND top100 | 0.629189 | 0.561817 | 0.714922 |

The analysis used 100,000 shared family draws, PCG64 seed 20260923,
and Bonferroni adjustment over six endpoints, without changing the frozen
protocol. F1 is the harmonic mean of macro precision and macro recall,
not mean family F1 or a pooled-pair statistic. Native half-counts and the
one-count prior are retained. Reconstructed and rounded native endpoints
agree within the prespecified tolerance.

For all-hit minus HMM, F1 was -0.061092 with adjusted interval
[-0.188023, 0.028017]; precision was -0.082110 [-0.181914, -0.001387];
recall was -0.029194 [-0.247769, 0.114154]. Top100 showed the same pattern:
F1 -0.059995 [-0.186513, 0.029034], precision -0.082123
[-0.181927, -0.001404], and recall -0.026345 [-0.245355, 0.114662].

Thus the adjusted precision differences favor initial HMM, but neither F1
nor recall difference excludes zero. For both controls, family-level F1
wins/ties/losses versus HMM were 5/4/9. These descriptive family counts do
not determine the aggregate statistic. This is development-exposed evidence
conditional on 18 families, not independent validation. Disjoint genes do
not establish family exchangeability. Equal search E-values do not establish
matched sensitivity or cost, and these controls are not standalone tools.

## Evidence And Verification

- [Source-bound count audit](qfo_sequence_swiss_counts_20260918.json),
  SHA-256 `3f04b05435f5967630e6e109e91b9e51d8618d940eb594e07681f5de84d0ae8c`.
- [Production intervals](qfo_sequence_swiss_bootstrap_20260918.json),
  SHA-256 `e4e4696e6f090d7d83b99e5190bd21606a980fac294c62cdf260582d822d89f3`.
- [Full interval table](qfo_sequence_swiss_bootstrap_20260918.md), including
  nominal intervals and wins/ties/losses for every metric.
- [Frozen protocol](QFO_SEQUENCE_UNCERTAINTY_PROTOCOL_20260918.md) and
  [admission inventory](qfo_sequence_swiss_admission_inventory_20260918.json).

All 80 focused exporter, count-audit, bootstrap and production-runner tests
passed. A separate NumPy calculation, using the retained confusion counts
without importing project statistic helpers, reproduced all three scalar
point estimates and all 12 nominal/adjusted interval pairs within 1e-12.
The production runner itself rehashed sources and reconstructed the count
audit before admitting intervals. It did not rerun upstream inference or
scoring. No uncertainty for other QfO endpoints or the secondary mean is
established by this analysis; full factorial intervals still await R-on cells.

Reproduce from the repository root with the recorded input paths intact:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
python benchmark_tools/run_qfo_sequence_uncertainty.py \
  --counts benchmark_tools/results/qfo_sequence_swiss_counts_20260918.json \
  --counts-sha256 3f04b05435f5967630e6e109e91b9e51d8618d940eb594e07681f5de84d0ae8c \
  --protocol benchmark_tools/results/QFO_SEQUENCE_UNCERTAINTY_PROTOCOL_20260918.md \
  --output /tmp/qfo_sequence_swiss_reproduction.json \
  --markdown /tmp/qfo_sequence_swiss_reproduction.md
```

Output paths must be fresh; underlying source files are required. This
command is source-bound reproduction, not a portable raw-data release.
