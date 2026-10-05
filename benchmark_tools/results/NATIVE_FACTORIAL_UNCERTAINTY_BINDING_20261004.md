# Native OrthoBench Results Linked To Retained Uncertainty

The four terminal-scored profiles-off cells have exactly the same 70
family records each as the retained eight-cell factorial. This includes
family identities, sizes, TP/FP/FN, splits and exact-family metadata, not
just equal aggregate F1/P/R. Therefore the retained deterministic paired
RefOG intervals can be reused for contrasts whose two cells are both
matched. No new bootstrap draws or independent accuracy are claimed.

The [binding helper](../bind_native_orthobench_uncertainty.py) checks pinned
snapshot/evidence, exactly replays the reporting snapshot through its
unchanged exporter, preserves failed-wrapper recovery, then requires full
native/cached family-record and point-estimate equality. Any mismatch
refuses reuse. The original factorial is fixed to SHA256
`6a0d588b5cb47c60fc6bc8bae8aa0c83e5f2aadb11de970919d8c6527c387141`.
The [actual binding](native_factorial_uncertainty_binding_20261004_v2.json)
retains all 12 planned contrasts and the 36-endpoint adjustment. Four are
matched; eight have unavailable native inputs and no imputed intervals.
Their retained cached analyses still exist separately.

| Factor | On - Off | F1 Difference (pp) | Nominal 95% CI | Adjusted Percentile CI | Family Wins/Ties/Losses |
| --- | --- | ---: | --- | --- | --- |
| Candidate expansion | P0C1R0 - P0C0R0 | -3.124 | [-7.568, 0.896] | [-10.659, 3.275] | 17/25/28 |
| Candidate expansion | P0C1R1 - P0C0R1 | 0.697 | [-1.396, 3.037] | [-2.459, 4.821] | 19/30/21 |
| Reconciliation | P0C0R1 - P0C0R0 | 2.942 | [0.841, 5.949] | [0.254, 8.144] | 11/59/0 |
| Reconciliation | P0C1R1 - P0C1R0 | 6.763 | [3.159, 11.285] | [1.376, 14.612] | 23/42/5 |

These reuse the original 20,000 shared paired draws, seed 20260918 and
alpha 0.05, with Bonferroni tails over all 36 planned metric endpoints.
Exact F1/P/R intervals are retained in JSON. Weighted F1 is recomputed
within each original draw; it is not a mean of per-family F1. Family
wins are descriptive, not another inferential endpoint.

This is development-exposed component evidence. Exchangeability of
RefOGs can fail through shared history or fused predictions; percentile
and multiplicity-adjusted coverage is approximate. Reconciliation changes
candidate co-membership to final root-HOG co-membership, not resolved-pair
truth. Profile-off retains sensitive initial HMM search. The intervals
neither isolate the whole HMM contribution nor establish superiority over
OrthoFinder. Failed22427 is still FAILED1:0 with scientific recovery, not a
clean-success timing. No timing comparison or correction is made here.

## Validation

The [independent stdlib readback](native_factorial_uncertainty_readback_20261004.json)
checks all 280 family records against the frozen source, independently
recomputes weighted F1/P/R with summation tolerance 1e-10 percentage points,
checks all four exact interval/metadata/count projections and eight null
native contrasts, preserves failed-wrapper status and verifies all 920
unchanged helpers. It does not re-review transitive raw artifacts or execute
another bootstrap. Source identity checks are bounded, not continuous
runtime closure.

All 68 joined tests pass, zero failures/errors/skips in 0.92s, including
16 new projection tests for changed family counts/metadata despite equal
aggregates, reordered records, missing/duplicate families, wrong cells,
changed multiplicity scope, absent contrasts and unavailable native input.
The [JUnit](native_factorial_uncertainty_binding_tests_20261004_v2.xml)
is retained. The first local binding preceded the final exporter-source
pin; v2 is a fresh non-overwriting binding with that check, not a changed
scientific analysis. Older local outputs are preserved.

## Reproduce

Run from the repository root in the existing analysis environment:

```bash
benchmarks/work/release_alert_refresh_20261001/venv/bin/python -B \
  -m benchmark_tools.bind_native_orthobench_uncertainty \
  --snapshot benchmark_tools/results/native_factorial_progress_20261004_v4/report.json \
  --snapshot-sha256 39369ecbc99f09dd2fe56a31a63bb0ffb4f52c3a40e20bbd9197af45fe824496 \
  --factorial benchmark_tools/results/orthobench_factorial_results_20260916.json \
  --output /tmp/native_orthobench_uncertainty_fresh.json
```

The output must not already exist. Original pinned local receipt paths
must be available; this is not a portable raw-data archive certification.
Remaining full-native cells, generalization, uncertainty and publication
distribution requirements remain active. Do not infer publication readiness.
