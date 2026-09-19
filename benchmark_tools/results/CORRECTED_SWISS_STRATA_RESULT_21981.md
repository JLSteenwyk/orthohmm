# Corrected SwissTrees Composition Strata

## Recovery And Provenance

Original job `21896` failed with exit 1:0 after seven seconds, before the
bootstrap call. Its driver independently reconstructed raw counts, then
incorrectly expected unseparated cell labels (`p1c0r0`) instead of the
canonical labels (`p1_c0_r0`). The synthetic test fixture repeated the same
mistake. The original executor and failed log are preserved; log SHA-256:
`a226f13e877b68c9c3afb48b4db51a597e9ff9fad442dfb0adb8e0cafcada1b4`.

Commit `6f81958c1d5fae5d39fa17beff0b332406e22147` imports the existing
hash-pinned canonical cell inventory and adds a regression test against the
upstream count-auditor contract. It changes no data, inference, reference,
stratum, bootstrap kernel, random seed, endpoint or multiplicity rule.
All 74 driver/kernel/descriptor/count-auditor tests passed. The unchanged
protocol SHA-256 remains
`2c87ba1df8ee39dfc325eefcdf374e713ba8fab1279a0330c3ea600ab637da54`.

Retry `21981` completed exit 0:0 in seven seconds with two CPUs/64 GiB on
`bizon`. The runtime log confirms Python 3.10.13 and NumPy 2.2.6. The frozen
executor is `benchmarks/work/qfo_corrected_strata_executor_6f81958`.
An initial submission was rejected for a scheduler dependency problem.
Before submitting without that dependency, accounting independently
confirmed both prerequisites (`21894` and `21736_0`) COMPLETED with exit
0:0. No incomplete prerequisite was bypassed; the driver still reconstructs
and checks the hash-bound evidence. Submission:

```bash
sbatch --parsable --dependency= \
  benchmark_tools/results/qfo_corrected_swiss_strata_retry_batch_20260919.sh \
  /mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_corrected_strata_executor_6f81958 \
  6f81958c1d5fae5d39fa17beff0b332406e22147
```

Driver SHA-256:
`93fdb5362ac387d4849214f6b99140c64dc24ac0b42b0375a25cd2da77343a50`.
Retry batch SHA-256:
`3eeaf61c78b773fb2afac873e97fcc6865b141402a230e89743465a79b71a3f0`.
The fresh output `benchmarks/work/qfo_corrected_swiss_primary_strata_v2.json`
is retained byte-for-byte as
[the result](qfo_corrected_swiss_primary_strata_21981.json), 278,164 bytes,
SHA-256 `93854aeddfd5243f407b0065b0a8bed74fbf3fc2a37775baca7c6301f34aa8ce`.

## Results

All 18 families are represented, nine in each frozen entropy bin and none
in the missing bin. The prespecified 100,000 paired family resamples use
PCG64 seed 20260924 and Bonferroni adjustment across 27 endpoints. F1 is
recomputed from macro precision/recall, not averaged per-family F1.

| Contrast | Entropy bin | F1 difference | Adjusted interval |
| --- | --- | ---: | --- |
| High sensitivity minus full OrthoFinder | Lower | -0.111298 | [-0.239322, 0.149647] |
| High sensitivity minus full OrthoFinder | Higher | -0.231373 | [-0.374440, -0.104983] |
| Phylogenetic minus full OrthoFinder | Lower | -0.028710 | [-0.154364, 0.147961] |
| Phylogenetic minus full OrthoFinder | Higher | -0.006102 | [-0.089129, 0.068670] |
| Phylogenetic minus high sensitivity | Lower | 0.082588 | [-0.020216, 0.181174] |
| Phylogenetic minus high sensitivity | Higher | 0.225271 | [0.067633, 0.396112] |

All adjusted intervals for phylogenetic OrthoHMM versus full OrthoFinder
include zero in both bins, across F1, precision and recall. High-sensitivity
OrthoHMM has lower precision in both bins with adjusted intervals excluding
zero; its higher-entropy F1 interval also excludes zero. Phylogenetic
OrthoHMM improves precision relative to high sensitivity in both bins and
F1 in the higher-entropy bin under this adjustment. All adjusted interaction
intervals and all adjusted recall intervals include zero (including a
recall lower bound exactly equal to zero). There are six nonzero adjusted
intervals among the 27 endpoints; all directions and null results are retained.

This does not demonstrate equivalence to OrthoFinder, independent
generalization, or a causal composition effect. The configuration contrast
changes candidate expansion as well as phylogeny. These are exploratory,
development-exposed conditional family intervals over small curated bins;
the adjustment does not cover previous development or other QfO metrics.

## Independent Check And Figure

[Numerical reproduction](qfo_corrected_swiss_strata_reproduction_21981.json)
independently recomputes the smoothed probabilities, harmonic macro F1,
paired family draws, all bin points/intervals, wins/ties/losses and all
interaction intervals. It uses family-wise weighted sums instead of the
production matrix product, with no production bootstrap/statistic helper
imports. Every endpoint agrees within absolute tolerance `1e-12`.
The check rehashes bound evidence before and after calculation. It reuses
NumPy's RNG/quantile implementation and admitted counts; it is not another
independent inference or raw-count audit. Its source is committed at
`71f7032`; all 64 reproduction/driver/kernel/plotter tests passed.

The [PNG](corrected_swiss_strata_figure_21981/corrected_swiss_strata.png),
[PDF](corrected_swiss_strata_figure_21981/corrected_swiss_strata.pdf), and
[SVG](corrected_swiss_strata_figure_21981/corrected_swiss_strata.svg) display
all 27 endpoints. The [TSV](corrected_swiss_strata_figure_21981/endpoints.tsv)
retains raw-scale values; the display alone multiplies by 100. The PNG was
visually inspected: no clipped intervals or overlapping labels. The figure
manifest binds input/source/output hashes. It does not claim publication
readiness. The all-method and secondary-stratum displays and complete
publication figure-bundle refresh remain unfinished.
