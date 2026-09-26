# Complete Descriptive SwissTrees Figures

Regenerated the identity, historical fragment and mapped-tree duplication
figures from the [complete admitted tables](SWISS_COMPLETE_STRATA_20260926.md).
The plotting entry points retain historical defaults; `--complete` selects
the explicitly pinned September 26 table hashes and requires all eight methods.

- [Identity PDF](swiss_descriptive_feature_figures_20260926/swiss_identity_descriptive.pdf): 48 displayed endpoints.
- [Fragment PDF](swiss_descriptive_feature_figures_20260926/swiss_fragment_descriptive.pdf): 96 displayed endpoints.
- [Duplication PDF](swiss_duplication_figure_20260926/swiss_duplication_descriptive.pdf): 48 displayed endpoints.

PNG, SVG, endpoint TSVs and source/output manifests accompany each figure.
All 192 displayed differences are available; empty bins remain in the source
tables. No missing value is converted to zero. OrthoMCL labels no longer say
unavailable in complete mode. Historical mode retains its missing-method
checks and labels. Twelve tests passed, including explicit-mode enforcement,
rejection of historical tables in complete mode, endpoint inventories,
delta/F1 validation and stale-label checks.

All three PNGs were visually inspected: labels, legends and notes fit,
OrthoMCL points render, and points are not clipped. No confidence intervals
or significance claims are added. Identity remains alignment-dependent,
unflagged is not proven complete, and duplication annotations are not an
evolutionary duplication rate or causal explanation.

```bash
python -m benchmark_tools.plot_swiss_descriptive_features --root . \
  --output /tmp/swiss-complete-feature-figures --complete
python -m benchmark_tools.plot_swiss_duplication_strata --root . \
  --output /tmp/swiss-complete-duplication-figure --complete
```

Output directories must be fresh. This reproduces plots from pinned tables,
not raw inference or feature admission. Manuscript links and preview updates
remain a subsequent integration step.
