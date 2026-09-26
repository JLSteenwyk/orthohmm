# Synthetic Family-Pair Variance Screen

Executed all eight cells and 8,000 replicates from the
[prespecified protocol](DYADIC_F1_VARIANCE_PROTOCOL_20260926.md), committed as
`d49c1a3` before execution. Numerical thread limits were one. No biological
benchmark outcomes entered simulation design or fitting; no benchmark
confidence intervals or scientific predictions changed.

## Coverage

Each cell has 1,000 replicates. All candidate variances were finite and
positive. Coverage counts retain every replicate, and Wilson Monte Carlo
intervals are available in the machine-readable result.

| Case | Families | Candidate 95% coverage | Naive coverage | Candidate variance / empirical variance |
| --- | ---: | ---: | ---: | ---: |
| Independent sparse | 64 | 0.937 | 0.943 | 0.951 |
| Independent sparse | 256 | 0.953 | 0.951 | 1.024 |
| Node sparse | 64 | 0.932 | 0.978 | 0.937 |
| Node sparse | 256 | 0.950 | 0.985 | 1.023 |
| Node dense | 64 | 0.934 | 0.981 | 0.985 |
| Node dense | 256 | 0.952 | 0.963 | 1.016 |
| Clade sparse, violated assumption | 64 | 0.791 | 0.881 | 0.440 |
| Clade sparse, violated assumption | 256 | 0.827 | 0.895 | 0.470 |

The six in-model cells meet the fixed necessary screen: no invalid variance
and coverage at least 0.925. This is not proof of nominal coverage, including
for all finite sample sizes. The naive independent-contribution calculation
is sometimes conservative here; dependence does not always make it too small.
The clade stress cases fail substantially for both calculations and were not
removed or used to select a correction.

## Independent Algebra Check

After the screen, a separate standard-library implementation derived the
exact variance of the *population-linearized* paired contrast under the
simulation model. It uses conditional Bernoulli/Poisson variance plus variance
of the conditional expectation over family effects. Groups of eight shared
effects multiply only the latter component by eight. This is a post-screen
diagnostic, not a changed acceptance rule or exact nonlinear F1 variance.

Across eight cells, empirical sample-F1 variance divided by this analytic
linearized variance ranges from 0.9727 to 1.0406. In the six in-model cells,
mean candidate variance / analytic variance ranges from 0.9599 to 1.0002.
In the two clade cells it is 0.4544 and 0.4749. Thus the stress-test failure
is consistent with omitted dependence, not merely a shifted simulated target.
This is evidence within the specified model, not proof of the cause of native
VGNC uncertainty. Twelve tests pass, including brute-force overlap algebra,
analytic derivatives, paired-method cancellation/symmetry, missing zero-dyad
rejection and an independent three-outcome Bernoulli variance enumeration.

## Retained Evidence

- [Screen](dyadic_f1_variance_screen_20260926.json), SHA-256
  `41088e1efe30d2661cfc8a7364c306ab318610c950d7809e8621242ec2eef66b`.
- [Analytic readback](dyadic_f1_variance_oracle_20260926.json), SHA-256
  `a03203fa56ea721a9800a722900fdbe905114efbd35356798b0ecb2a7ca01b49`.
- Sources: `validate_dyadic_f1_uncertainty.py` and
  `audit_dyadic_f1_screen.py`; each result records its source checksum.

Reproduce the simulation with the protocol command and a fresh output path.
The independent diagnostic is:

```sh
/usr/bin/python3 -S -m benchmark_tools.audit_dyadic_f1_screen \
  --source benchmark_tools/results/dyadic_f1_variance_screen_20260926.json \
  --output NEW_ORACLE.json
python -m pytest -q tests/unit/test_validate_dyadic_f1_uncertainty.py \
  tests/unit/test_audit_dyadic_f1_screen.py
```

## Decision and Next Requirement

The candidate merits further evaluation but is **not admitted for native
benchmark intervals**. The [underlying dyadic reference](https://arxiv.org/abs/1312.3398)
requires a dependence model that cannot be inferred from disjoint reference
labels alone. Singleton contributions, ratio behavior and the synthetic
extension are addressed only in this limited screen.

Before any VGNC application, map corrected native TP/FP/FN rows to fixed
reference blocks for all retained methods; preserve every denominator and
zero opportunity; assess family-size imbalance, small-denominator behavior,
and shared species/annotation/prediction dependencies. Define the biological
sampling target and any fixed higher-level sampling units independently of
method errors. Outcome-dependent graph components are not a solution. If no
defensible sampling interpretation is available, retain VGNC point estimates
and explicit uncertainty limitations rather than publish these synthetic
intervals as validated biological inference. GO/EC, sampled FAS and the
project-defined secondary mean need separate treatment.
