# Conditional Paired Poisson F1 Validation

All 12 prespecified conditional-model cells passed the numerical coverage check.
Both non-Poisson common-shock controls had coverage zero. The construction is
therefore a restricted mathematical ingredient, not an admitted native benchmark
uncertainty method or evidence of improved OrthoHMM accuracy.

## Evidence And Target

The [frozen protocol](PAIRED_POISSON_F1_PROTOCOL_20261010.md) specifies fixed,
perfect-recall truth mass T and independent whole-panel Poisson counts for
shared and method-only false positives. The target is F1 evaluated at the
expected count vectors, differenced between methods:

```text
delta = 2T/(2T+c+a) - 2T/(2T+c+b).
```

This is not the expectation of sample F1 or inference from a fixed native
reference table to unseen families. It does not assume that native gene pairs
or reference families are IID.

[Garwood limits](https://jhanley.biostat.mcgill.ca/bios601/Intensity-Rate/Garwood1936.pdf),
in the [PDG chi-square form](https://pdg.lbl.gov/2021/reviews/rpp2021-rev-statistics.pdf),
provide the component Poisson mean intervals. The project's union-bound and
sharp-rectangle projection argument is in the protocol. In ideal arithmetic,
valid Poisson marginals suffice for its at-least-95% coverage argument;
independence is additionally assumed for the joint numerical enumeration.

## Frozen Panel Results

The following table is generated from the retained
[validation JSON](paired_poisson_f1_validation_20261010_v1.json).
Covered mass is unnormalized probability mass summed over the truncated
enumeration, not a Monte Carlo coverage estimate. Mean width is likewise the
probability-weighted sum on that enumeration, not a native benchmark interval.

| Truth Mass | Means (Shared, Left Only, Right Only) | F1 Difference Target | Covered Enumerated Mass | Mean Interval Width |
| --- | --- | --- | --- | --- |
| 512 | (0, 0, 0) | 0.0000000000e+0 | 1.000000000000 | 9.3070566686e-3 |
| 512 | (2, 0, 0) | 0.0000000000e+0 | 1.000000000000 | 9.3023753640e-3 |
| 512 | (1, 1, 0) | -9.7370798269e-4 | 0.999999999936 | 1.1111982482e-2 |
| 512 | (1, 0, 1) | 9.7370798269e-4 | 0.999999999936 | 1.1111982482e-2 |
| 512 | (0.5, 0.3, 0.7) | 3.8986319152e-4 | 0.999999999988 | 1.1182854708e-2 |
| 512 | (10, 20, 12) | -7.4304933224e-3 | 0.999583309608 | 3.7560535431e-2 |
| 23934 | (0, 0, 0) | 0.0000000000e+0 | 1.000000000000 | 2.0000889832e-4 |
| 23934 | (2, 0, 0) | 0.0000000000e+0 | 1.000000000000 | 2.0000673886e-4 |
| 23934 | (1, 1, 0) | -2.0889473776e-5 | 0.999999999936 | 2.3916129465e-4 |
| 23934 | (1, 0, 1) | 2.0889473776e-5 | 0.999999999936 | 2.3916129465e-4 |
| 23934 | (0.5, 0.3, 0.7) | 8.3559640659e-6 | 0.999999999988 | 2.4064006022e-4 |
| 23934 | (10, 20, 12) | -1.6694486015e-4 | 0.999559608605 | 8.3944324154e-4 |

The minimum covered mass is 0.999559608605.
The maximum omitted probability mass is
5.999645225074e-13.
The intervals are conservative in this panel. Decimal display rounds the
numbers; tail-only coverage bounds do not certify floating-point rounding.

## Model Violation Controls

For each truth mass, the control sets shared and right-only counts to zero.
The left-only count is zero with probability 0.99 and 1,000 with probability
0.01, giving mean 10 but a non-Poisson marginal. Neither outcome's interval
covers the mean-count F1 target: coverage is exactly zero for both controls.
These failures remain in the JSON and are not covered Poisson regimes.

Imperfect recall, unequal-size latent-node laws, shared-clade dependence and
native VGNC are outside the construction. Earlier failed dyadic screens remain
unchanged; passing this panel does not replace them or justify biological
resampling units. Native VGNC, TreeFam, GO, EC, FAS and secondary-mean paired
intervals remain unresolved. No native endpoint or denominator was substituted.

## Execution And Reproduction

Tested source and protocol were committed and pushed as `e7dbf881` before
the actual enumeration. The [execution record](paired_poisson_f1_execution_20261010_v1.json)
retains the precise commands, environments and terminal results. Enumeration
and independent readback both exited zero. The latter independently inverted
Poisson tails and evaluated all eight rectangle corners, without importing the
producer; maximum discrepancy was 1.4432899320127035e-15.
The focused and adjacent suites passed all 60 tests in 2.52 seconds.

The retained runtime is Python 3.12.3, NumPy 2.2.6 and SciPy 1.15.3 with one
numeric thread. The command refuses occupied output paths. Reuse the retained
result rather than repeating this completed panel; any deliberate reproduction
must use a fresh output path and preserve the explicit frozen protocol digest.

Result SHA256:
`314ddd68b690fb8a8ef210f49fa5be85e3a1a6a17b8f6362cefbf9d091ebee07`.
The execution record and result bind source and protocol bytes. This integration
changes reporting only: no new inference, scoring, timing, random draws,
parameter tuning, independent biological confirmation or publication admission.

