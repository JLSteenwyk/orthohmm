# Corrected QfO Factorial Figure

Rendered the independently validated corrected-release SwissTrees factorial
without resampling, selecting endpoints or changing the underlying results.
The four-panel figure contains all eight cell scores and all 42 contrast
endpoints (14 contrasts each for F1, precision and recall). Both nominal
and adjusted intervals are shown in percentage-point units.

- [PNG](qfo_corrected_factorial_figures_20260919/qfo_factorial_swiss.png)
- [PDF](qfo_corrected_factorial_figures_20260919/qfo_factorial_swiss.pdf)
- [SVG](qfo_corrected_factorial_figures_20260919/qfo_factorial_swiss.svg)
- [Manifest](qfo_corrected_factorial_figures_20260919/manifest.json)

The shared plotter requires explicit `--input-release corrected`, the
corrected result status, exact input-release identity and frozen corrected
protocol hash. Historical and corrected results cannot be silently swapped.
The caption derives the number of adjusted F1 intervals excluding zero
from the supplied results (4/14), not a hard-coded scientific conclusion.
The original default remains historical and existing figures are unchanged.

Twenty-three focused tests cover both releases, all 42 plotted point
estimates and interval coordinates for the corrected release, score units,
cell/contrast completeness, protocol mutation and release confusion. Figure
text bounds pass; visual inspection found no clipped or overlapping labels.
The PNG is 3420 by 1620 pixels with 11.99% nonwhite pixels. All manifest
source, helper and output hashes were rechecked after rendering.

Source result SHA-256:
`f777dead1294b3810c0479877aade7fb4411c0544696f671226ff198432e9211`.
Figure manifest SHA-256:
`133fc7bb1d2404b7c9e5d932486c6134ff1f2fa44e1880dd839e7fa7fe3e4134`.

The figure is linked in the manuscript and preserved in a
[separate corrected-release evidence supplement](PUBLICATION_CORRECTED_FIGURE_BUNDLE_20260919.md).
The earlier bundle is unchanged; a complete publication archive remains unfinished.
This is development-exposed evidence from 18 families, not independent
confirmation, selection-adjusted inference or superiority over OrthoFinder.
P-off retains initial HMM search, and R changes prediction semantics.

## Reproduction

```bash
python -m pytest -q tests/unit/test_plot_qfo_factorial.py tests/unit/test_plot_corrected_qfo_factorial.py
python -m benchmark_tools.plot_qfo_factorial \
  --results benchmark_tools/results/qfo_corrected_factorial_complete_20260919/swiss_bootstrap.json \
  --sha256 f777dead1294b3810c0479877aade7fb4411c0544696f671226ff198432e9211 \
  --input-release corrected --output /new/path/corrected_factorial_figure
```
