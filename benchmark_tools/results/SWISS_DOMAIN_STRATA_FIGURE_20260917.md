# SwissTrees Domain-Stratified Figure

[PNG](figures/swiss_domain_strata_20260917/swiss_domain_strata.png),
[PDF](figures/swiss_domain_strata_20260917/swiss_domain_strata.pdf),
[SVG](figures/swiss_domain_strata_20260917/swiss_domain_strata.svg).

## Caption

Paired method differences in SwissTrees F1, precision and recall, stratified
by family median distinct Pfam types. Lower and higher bins contain 12 and
six families, respectively. Each comparison also shows the higher-minus-lower
interaction. Points represent observed differences in percentage points;
thick lines indicate nominal 95% percentile intervals and thin lines indicate
Bonferroni-adjusted intervals across all 27 prespecified endpoints. The frozen
analysis used 100,000 within-bin family resamples, paired across methods,
PCG64 seed 20260921, recomputing macro precision/recall and harmonic F1.

Both OrthoHMM configurations have negative adjusted F1 differences versus
full OrthoFinder in both bins. All nine adjusted interaction intervals include
zero: a positive result in one bin but not the other does not establish a
difference between bins. This is retrospective, development-exposed evidence
with possible family dependence, not a causal domain-effect analysis. The
OrthoHMM configurations differ beyond reconciliation. Repeat-type bins are
descriptive only and are tabulated separately, not plotted here. Intervals
do not cover other QfO endpoints or the secondary six-metric mean.

## Provenance and Validation

The plotter requires the exact frozen result SHA-256
`229269f3db0bfb0e99893de00a1c2467285a07d7ca148c4854b7661ca13a5891`.
No statistics were recalculated or thresholds changed. The figure manifest
records input, plotter and all output hashes plus NumPy/Matplotlib versions.
Generation used the patched, pinned Swiss-analysis environment.

Seven new plotting tests check all 27 points and 54 intervals against source
values, rendered text bounds, and rejection of missing comparisons, nonfinite
values, reversed intervals, axis overflow, changed bin sizes and changed row
order. Together with the existing analysis and comparator-plot tests,
17 tests pass. The generated PNG was inspected for legibility and overlap.

From the repository root, use a new output directory:

```bash
benchmarks/work/swiss_analysis_env_20260917/bin/python benchmark_tools/plot_swiss_domain_strata.py --results benchmark_tools/results/swiss_domain_strata_results_20260917.json --output /tmp/swiss-domain-figure
python -m pytest -q tests/unit/test_plot_swiss_domain_strata.py tests/unit/test_analyze_swiss_domain_strata.py tests/unit/test_plot_qfo_swiss_comparators.py
```

This figure completes a presentation artifact, not the remaining empirical
QfO factorial, resource admission or publication-readiness requirements.
