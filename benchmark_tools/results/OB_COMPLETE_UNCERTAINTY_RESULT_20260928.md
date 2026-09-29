# Eight-Method OrthoBench Paired Uncertainty

Protocol and runner were pushed as `6d25b18f` before execution. The analysis
retains all 70 RefOGs, eight methods and full-reference weighted counts.
100,000 paired multinomial draws use PCG64 seed 20260928; all seven methods
are contrasted with full OrthoFinder 3.1.5 for F1, precision and recall.
Bonferroni percentile intervals cover 21 planned endpoints, with approximately
119 draws per adjusted tail. They are approximate, conditional intervals.

[Complete 21-endpoint table](ob_complete_uncertainty_check_20260928_v2/TABLE.md)
and [full result](ob_complete_uncertainty_20260928.json) retain every estimate,
nominal/adjusted interval and descriptive per-family F1 win/tie/loss count.

| Method | F1 difference vs full OF (pp) | Adjusted interval |
|---|---:|---|
| OrthoHMM high sensitivity | -2.377 | [-11.045, 6.746] |
| OrthoHMM phylogenetic | +1.370 | [-7.281, 12.166] |
| OrthoFinder sequence-only | -14.031 | [-26.657, -1.859] |
| SonicParanoid | -25.979 | [-52.286, 5.451] |
| ProteinOrtho | -27.679 | [-39.688, -12.857] |
| FastOMA | -41.830 | [-54.989, -25.918] |
| OrthoMCL | -17.671 | [-32.754, -1.040] |

Neither OrthoHMM F1 interval excludes zero. Phylogenetic OrthoHMM retains
an adjusted precision advantage (+15.705 points, [1.175, 30.634]) and recall
deficit (-13.151, [-26.947, -0.609]). High-sensitivity precision no longer
excludes zero under this broader 21-endpoint correction; its recall deficit
remains. These observations neither prove equivalence nor population-wide
superiority. Preserve the original six-endpoint comparison as historical;
do not substitute its narrower intervals into this expanded panel.

An [independent accumulation check](ob_complete_uncertainty_check_20260928_v2/crosscheck.json)
reproduced all point estimates and 21 nominal/adjusted intervals, with maximum
error 2.14e-14 percentage points. It sums each family's contribution explicitly
and uses exact rational observed counts instead of the primary matrix product.
RNG and quantile implementations remain shared NumPy code. The v2 check also
explicitly verifies reported point differences; the first check is retained
locally. This is arithmetic agreement, not validated interval coverage.
Twenty-four focused tests pass. Native inference/scoring were not rerun.

Family exchangeability is an assumption, not established biological independence.
Shared histories, shared genes and errors spanning families may violate it.
The benchmark and method selection are development-exposed. No tuning,
pair-IID intervals, bootstrap p-values or new independent confirmation is
introduced. Remaining QfO uncertainty and publication gates stay open.

```sh
python -B -m benchmark_tools.bootstrap_ob_complete --base benchmark_tools/results --output /tmp/ob-complete-uncertainty-new.json --protocol-sha256 a75847ff383b87a288e18953e8fce7858a4f939f175c428a720d5f5e3118cea3
python -B -m benchmark_tools.check_ob_complete_uncertainty --base benchmark_tools/results --output /tmp/ob-complete-uncertainty-check-new
```

Use fresh output paths. The second command checks the retained result in the
base directory, not the first command's `/tmp` output. Historical input records
refer to local files; this is not a standalone portable inference recipe.
