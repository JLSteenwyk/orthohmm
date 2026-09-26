# Synthetic Screen for Family-Pair F1 Uncertainty

This protocol precedes the 1,000-replicate execution. No benchmark score is
used to select settings. The candidate is not yet authorized for native VGNC,
GO/EC, FAS or the six-endpoint mean. Earlier audits established that many VGNC
false positives join different reference blocks; assigning them arbitrarily
to one family does not preserve the required dependence structure.

## Candidate and Assumptions

Use paired TP/FP/FN arrays indexed by a single family (within-family rows)
or an unordered pair of families (cross-family rows). All zero dyads remain
in the opportunity grid. Recompute each method's aggregate F1 from its total
counts. Linearize the paired difference with the analytic F1 gradients, using
the same simulated observations for both methods. Center within-family and
cross-family contributions separately, since their distributions differ.

If `d_i` is the centered within-family contribution and `e_ij` the centered
cross-family contribution, evaluate

```
V = sum_i (d_i + sum_{j != i} e_ij)^2 - sum_{i < j} e_ij^2
```

This is the algebraic sum of products for observations sharing at least one
family, counting each edge's self-product once. A brute-force unit test checks
that identity, including diagonals. It does not prove statistical validity.
Use an uncorrected normal 95% Wald interval for the paired F1 difference;
retain zero/negative/nonfinite variance as failure, not a clipped interval.
Compare against the naive sum of independent contribution squares.

The motivation is the shared-endpoint sandwich construction of
[Aronow, Samii and Assenova (2015)](https://arxiv.org/abs/1312.3398).
Their basic setting samples independent units, excludes self-dyads and
requires zero covariance for disjoint dyads. Our singleton contributions,
stratum centering, sparse-count models and nonlinear ratio require separate
validation; their theorem is not cited as proof of this extension. Biological
families, shared species histories and joint prediction fitting may violate
the independence assumptions even if simulations pass.

## Frozen Simulation Design

For each of 64 and 256 families, execute 1,000 independent replicates in each
of four cases. PCG64 uses `SeedSequence([20260927, case_index, family_count])`.
Each family has 20 true pairs. Methods A/B have mean recall 0.75/0.80 and
share Bernoulli uniforms. In heterogeneous cases, a symmetric +/-1 family
effect changes recall by 0.15/0.10. False positives occur between families;
shared Poisson draws give means 8/5 times a positive endpoint-effect factor,
divided by family count in sparse cases or by 32 in the dense case.

1. Independent sparse: no latent family heterogeneity.
2. Node sparse: independent family effects and sparse false positives.
3. Node dense: independent effects and dense false positives.
4. Clade sparse: groups of eight families share one effect. This deliberately
   violates independence of disjoint dyads; it is a stress test, not a gate.

The known target is the difference of F1 ratios evaluated at exact expected
count totals for the chosen finite family count. It is not the expectation of
sample F1. Report empirical bias/variance, mean candidate/naive variances,
interval coverage and Wilson Monte Carlo intervals, and invalid variances.
No cell, replicate or negative variance may be omitted or retried.

A necessary exploratory screen is zero invalid variances and at least 0.925
coverage in every in-model cell. This fixed tolerance is not a claim of exact
95% coverage. Failure prevents application to benchmark intervals without a
new justified protocol. Passing also does not establish general validity:
native denominator/eligibility mapping, reference exchangeability, dependence
outside shared families, sparsity, unequal family sizes, small denominators,
and multiplicity would still need review. No correction is chosen after
inspecting coverage, and no default or scientific prediction changes.

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python -m benchmark_tools.validate_dyadic_f1_uncertainty --replicates 1000 \
  --output benchmark_tools/results/dyadic_f1_variance_screen_20260926.json
python -m pytest -q tests/unit/test_validate_dyadic_f1_uncertainty.py
```
