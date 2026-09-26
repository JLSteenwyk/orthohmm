# Reference-Size and Rare-Error Variance Screen

Freeze this extension before its coverage run. It does not change the earlier
screen or authorize intervals for native benchmarks. The preceding block
mapping found 16,844 reference blocks, 14,075 with only one asserted pair;
the previous fixed-20-pair simulation did not represent that size distribution.

Read only the reference-pair-count histogram from the pinned mapping report
`01a51a62a536bcf759f2a0ef83c56196a92022b0bb8290448eacbfefa47d6005`.
Sample family true-pair counts independently from this empirical distribution.
Do not fit any parameters to tool-specific TP/FP/FN counts. Numerical cases
below are generic stress regimes, not fitted models of named tools.

## Frozen Design

Run 1,000 replicates for each combination of 256 or 16,844 families and:

1. Unequal regular: independent families, mean recalls .75/.80, expected
   total FP 4(n-1)/2.5(n-1), paired Bernoulli and Poisson draws.
2. Unequal node: symmetric independent family effects alter recalls by
   .15/.10 and FP intensity by the same additive endpoint factor as the
   earlier protocol. Conditional intensities are sampled by Poisson thinning.
3. Rare perfect: both recalls equal one, expected total FP two/one, irrespective
   of family count. The methods share the Poisson component with mean one.
   This deliberately probes boundary behavior with a non-growing error count.

PCG64 uses `SeedSequence([20260928, case_index, family_count])`. Every family
count is positive; the analytic target is the difference of F1 ratios at
expected totals, with expectation taken over sampled family sizes as well
as outcomes. The earlier normal-Wald estimator, threshold .925, and policy
of zero/negative/nonfinite variance counting as noncoverage are unchanged.
All six cells must pass for this broader screen to pass. A failure means no
blanket application; it does not prove every interior parameter setting fails.

## Sparse Algebra and Scope

For n families, M=n(n-1)/2 off-diagonal cells exist, including zeros. Center
diagonal contributions d by their own mean. For stored raw nonzero edge
contributions e, let m=sum(e)/M. Then compute the earlier variance without
materializing the full grid:

```
incident_i = centered_d_i + sum_j e_ij - (n-1)*m
V = sum_i incident_i^2 - (sum_stored_edges e_ij^2 - M*m^2)
```

Repeated simulated FP events are aggregated into unique unordered dyads
before variance calculation. Tests compare this expression with the previous
full-grid implementation for every case, including empty edge tables.
Sparse storage is an algebraic optimization, not a statistical correction.

No clade-dependence repair is introduced; the previous failure remains.
No multiple correction strategies will be selected after seeing coverage.
Retain bias, empirical variance, mean candidate variance, coverage with
Monte Carlo intervals and invalid counts. No benchmark CIs, native prediction
changes, DGX work or controlled timing claim is included.

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python -m benchmark_tools.validate_sparse_dyadic_f1 --repo . --replicates 1000 \
  --output benchmark_tools/results/sparse_dyadic_f1_screen_20260926.json
python -m pytest -q tests/unit/test_validate_sparse_dyadic_f1.py
```
