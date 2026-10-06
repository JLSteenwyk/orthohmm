# Native QfO Figure Protocol

`benchmark_tools.plot_native_qfo_scientific_scores` consumes the actual combined
scientific snapshot and its guarded SwissTrees interval binding. Directly replay
the existing report/binding checks, without raw count, inference, scoring or
bootstrap repetition. Refuse a changed snapshot, unadmitted plotted cell,
relabeled failed timing, missing endpoint, nonfinite value, changed contrast or
invalid interval. Require fresh outputs.

Plot P0/C0/R0 versus P0/C0/R1 only. Initial HMM search remains on and this is
not a selected-default or whole-tool comparison. Separate three orthology F1
endpoints from GO/EC/FAS similarity, show SwissTrees precision/recall and the
three R-effect intervals. Retain nominal and all-42-endpoint adjusted intervals,
including the adjusted F1 zero crossing. Do not draw native SEM error bars,
substitute unavailable cells or add a custom aggregate to the primary figure.

Export full-precision score and interval TSVs plus PNG/PDF/SVG and an input,
source, package-version and output-hash manifest. The initial manifest explicitly
leaves visual review incomplete; inspect rendered assets separately. This is
not a new score admission, native uncertainty calculation, timing comparison,
hermetic environment certification or completed publication package.

129 joined tests pass in4.95s, including25 new plot cases. Test fixture rendering
checks2500x1600 nonblank pixels, both method colors, SVG labels, endpoint
exports, no overwrite and scope/failure handling. The fixture stubs direct
binding replay; only the subsequent actual invocation checks real provenance.
Matplotlib3.10.8/Pillow12.3.0 already exist in the retained Python3.12 environment;
no package installation or scientific runtime change.

```bash
benchmarks/work/release_alert_refresh_20261001/venv/bin/python -B \
  -m benchmark_tools.plot_native_qfo_scientific_scores \
  --snapshot benchmark_tools/results/native_qfo_scientific_scores_20261006_v1/report.json \
  --snapshot-sha256 6b2f735ea8a44f72715e11a7f015e4328576baa74c544889cbed0cdfe70ea07b \
  --swiss-binding benchmark_tools/results/native_qfo_swiss_uncertainty_binding_22449_20261006.json \
  --swiss-binding-sha256 85673a114b7f7c6e05de5da189c8cbe8d77ec5d1378604dbae7e5902dfd996cc \
  --output benchmark_tools/results/native_qfo_p0c0_figure_20261006_v1
```

At preparation no actual figure is generated. Original native22444 remains
live and22445 dependency-pending. Preserve frozen scientific settings, prior
archive/PDF bytes, all other uncertainty requirements and shared-host caveats.
