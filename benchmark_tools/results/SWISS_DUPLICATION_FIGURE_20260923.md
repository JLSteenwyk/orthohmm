# Duplication-Annotation Figure

[Figure](swiss_duplication_figure_v2_20260923/swiss_duplication_descriptive.png)
and [machine-readable provenance](swiss_duplication_figure_v2_20260923/manifest.json)
derive directly from the pinned 32-row descriptive score table. The existing
validator checks complete method/bin coverage, availability, family counts,
finite values, harmonic F1 and differences against full OrthoFinder. The
older identity/fragment renderer and its pinned source remain unchanged.

All eight methods appear in each of three metric panels. Two nine-family
strata give 48 endpoint records, six unavailable for OrthoMCL. Shared axis
limits are -45 to +20 percentage points. No intervals or significance marks
are added. The empty missing-feature bin and full-panel result remain in
the source table; the figure does not invent points for them.

The initial plot limit of +15 rejected an out-of-range recall endpoint
before writing an image. Expanded to +20 and rendered into a fresh directory;
all endpoints pass bounds checks. The 2700-by-1620 PNG was inspected visually:
labels and footnotes are readable, points are visible, and no text overlaps
the plot or adjacent elements. PDF and native Matplotlib SVG are also retained;
only the PNG received visual inspection. Native SVG whitespace is preserved.
Eleven plot/validator tests pass, including nonblank pixel checks and changed
source rejection.

Reproduce with a fresh output directory:

```bash
python -m benchmark_tools.plot_swiss_duplication_strata --root . --output /tmp/swiss-duplication-figure
python -m pytest -q tests/test_plot_swiss_duplication_strata.py tests/unit/test_plot_swiss_descriptive_features.py
```

This is visualization of existing descriptive results, not new scientific
inference, a controlled resource benchmark or completion of the publication.
