# Corrected QfO Threshold Neighborhood

All four score jobs 22043_0..3 and independent score admissions 22047_0..3
completed with exit 0:0. The frozen analysis reconstructed native SwissTrees
family counts from raw predictions, checked reference identities and native
score arithmetic, and freshly verified the linked input identities.

## Results

SwissTrees F1 is the harmonic mean of macro precision and macro recall, with
the native per-family prior. Values below are percentages; differences and
intervals are percentage points relative to the frozen full-pipeline control.

| Arm | Changed parameter | F1 (%) | Difference | Adjusted interval | Family wins/ties/losses |
| --- | --- | ---: | ---: | --- | --- |
| Control | None | 83.3513 | 0 | Reference | Reference |
| norm_low | min_norm=0.024 | 83.3513 | 0 | [0, 0] | 0/18/0 |
| norm_high | min_norm=0.036 | 83.3513 | 0 | [0, 0] | 0/18/0 |
| margin_low | min_margin=1.2 | 82.8391 | -0.5122 | [-4.1173, 1.3523] | 1/16/1 |
| margin_high | min_margin=1.8 | 83.5648 | +0.2134 | [-0.9285, 1.9571] | 1/16/1 |
| cpm_low | cpm_resolution=0.08 | NA | NA | NA | NA |
| cpm_high | cpm_resolution=0.12 | NA | NA | NA | NA |

Control settings are min_norm=0.03, min_margin=1.5 and cpm_resolution=0.1.
Both normalization variants have identical family statistics to the control
on these 18 families; zero-width empirical intervals do not establish
equivalence on unobserved families or identical whole-proteome predictions.
Neither margin contrast excludes zero, including for precision and recall.
No frozen defaults were changed based on these development-exposed results.

Intervals use the prespecified 100,000 paired family bootstrap replicates,
PCG64 seed 20260925, linear quantiles, and Bonferroni correction across all
18 planned endpoints (six contrasts times F1, precision and recall). Missing
CPM contrasts remain in the correction. Low-CPM phylogeny is running;
high-CPM native replay failed with SIGSEGV. No missing scores are imputed.
The panel remains incomplete, and this is SwissTrees uncertainty, not a
joint uncertainty analysis of all QfO endpoints.

## Evidence and Reproduction

- [Admission inventory](qfo_parameter_uncertainty_threshold_inventory_20260923.json),
  SHA256 `829cdc5fd804bbd6e9c5e6b7166c51973c30625dcb5df491aac8ad7963ad8fe3`.
- [Count audit and bootstrap](qfo_parameter_uncertainty_threshold_20260923.json),
  SHA256 `0f5dafd4da19a4ed412f33af3ca9122dde37750f82f960712534f917192e292f`.
- [Independent arithmetic reproduction](qfo_parameter_uncertainty_threshold_reproduction_20260923.json):
  all 12 estimated endpoints agree within absolute tolerance 1e-12.

From the repository root, with retained source/input paths available, use
fresh output paths (existing reports are deliberately not overwritten):

```bash
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /home/bizon/anaconda3/bin/python benchmark_tools/run_qfo_parameter_uncertainty.py \
  --inventory benchmark_tools/results/qfo_parameter_uncertainty_threshold_inventory_20260923.json \
  --inventory-sha256 829cdc5fd804bbd6e9c5e6b7166c51973c30625dcb5df491aac8ad7963ad8fe3 \
  --baseline benchmark_tools/results/qfo_swiss_counts_20260917.json \
  --plan benchmark_tools/results/qfo_parameter_neighborhood_plan_20260919.json \
  --protocol benchmark_tools/results/QFO_PARAMETER_NEIGHBORHOOD_PROTOCOL_20260919.md \
  --output /tmp/qfo_threshold_uncertainty_fresh.json
```

The numerical reproducer is `benchmark_tools/reproduce_qfo_parameter_uncertainty.py`
with `--results`, `--results-sha256` and a fresh `--output`. It reconstructs
statistics and bootstrap differences independently from admitted counts,
but uses the same NumPy generator and quantile implementation; it does not
rerun upstream inference or replace the raw-count audit. All 117 focused
parameter-audit/bootstrap/reproduction tests pass. Controlled runtime
comparison remains deferred and is not inferred from these shared-host jobs.
