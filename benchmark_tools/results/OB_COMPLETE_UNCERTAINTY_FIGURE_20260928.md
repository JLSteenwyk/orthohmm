# Complete OrthoBench Comparison Figure

[PDF](figures_ob_complete_uncertainty_20260928/ob_complete_uncertainty.pdf),
[PNG](figures_ob_complete_uncertainty_20260928/ob_complete_uncertainty.png),
[SVG](figures_ob_complete_uncertainty_20260928/ob_complete_uncertainty.svg) and
[plotted values](figures_ob_complete_uncertainty_20260928/endpoints.tsv).

Seven candidate methods are compared with full OrthoFinder 3.1.5 for F1,
precision and recall. Dots are differences in percentage points; thick lines
are nominal 95% paired percentile intervals and thin lines retain the complete
21-endpoint Bonferroni adjustment. All panels share a fixed -80 to +60 scale
and method order. Full OrthoFinder is the zero reference, not an eighth
self-comparison. The sequence-only checkpoint remains labelled as diagnostic.

The figure imports the existing complete analysis; it does not rerun inference,
bootstrap sampling or endpoint selection. Neither OrthoHMM F1 interval excludes
zero. Phylogenetic OrthoHMM has an adjusted precision advantage and recall
deficit. Family exchangeability and development exposure limit interpretation;
zero inclusion does not prove equivalence. No overall superiority is claimed.

Ten plotter tests pass, covering the complete panel, retained percentage-point
units, missing methods, wrong baseline/correction, invalid intervals, nonfinite
values, inconsistent differences, checksum failure and overwrite refusal.
[Independent export check](ob_complete_uncertainty_figure_check_20260928.json)
compares all 105 exported numeric values directly to the retained JSON and
finds exact agreement; the one-page PDF has no out-of-bounds words. Both the
generated PNG and PDF raster were visually inspected with no clipping or
incoherent overlaps. This verifies presentation, not interval coverage or raw
prediction correctness. The automated render manifest does not claim visual
review or publication readiness.

```sh
python -B -m benchmark_tools.plot_ob_complete_uncertainty \
  --results benchmark_tools/results/ob_complete_uncertainty_20260928.json \
  --sha256 f472916ccafa2392b7ca27f39448195b7c3e6b9464772e0ca6775092908f6b75 \
  --output /tmp/ob-complete-uncertainty-figure-new
```

Use a fresh directory. Source and input checksums are retained in the figure
manifest. Prior six-endpoint figures remain historical and are not substituted
for this complete panel. Scientific source, timing workflow and defaults are
unchanged; remaining publication requirements remain open.
