# Complete Corrected-QfO Endpoint Figure

[PNG](figures_corrected_qfo_endpoints_20260926/corrected_qfo_endpoints.png),
[PDF](figures_corrected_qfo_endpoints_20260926/corrected_qfo_endpoints.pdf),
[SVG](figures_corrected_qfo_endpoints_20260926/corrected_qfo_endpoints.svg), and
[manifest](figures_corrected_qfo_endpoints_20260926/manifest.json).

All eight admitted methods and six endpoints are shown from the complete
corrected comparison manifest, pinned by SHA256
`042aa221d01554114a1b0b413ca3e2ff56dad5f8c06362c69a1bf49f2de642fc`.
The figure preserves the established method colors/markers and supplied-tree
FastOMA label. Historical original-release figures are unchanged.

VGNC, SwissTrees and TreeFam-A use precision-recall axes; the plotter verifies
that their harmonic means match retained F1 values. GO/EC use assessed relation
counts and mean functional similarity. FAS uses the reported eligible relation
count, not the size of its unseeded sample. Prediction semantics are retained
per method in the figure manifest. Forty-eight points are recorded; near-equal
points may overlap. No jitter or invented uncertainty was introduced.

Eight focused tests pass for inventory, admission, finite values, statistic
identity, F1 consistency, real rendering and overwrite refusal. The PNG was
visually inspected for clipping and layout. This is retained-score visualization,
not raw-score revalidation, uncertainty estimation, or a superiority claim.

```bash
python -m benchmark_tools.plot_corrected_qfo_endpoints \
  --source benchmark_tools/results/qfo_corrected_comparison_20260926_v7/manifest.json \
  --output /tmp/corrected-qfo-endpoint-figure-new
python -m pytest -q tests/unit/test_plot_corrected_qfo_endpoints.py
```

The working manuscript still links the original-release figure in its historical
section. Integration of this additional corrected figure into the corrected
section and a refreshed manuscript preview remains next. This artifact is not
a full release or publication-readiness declaration.
