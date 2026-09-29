# All-Method OrthoBench Figure

[PDF](ob_complete_strata_figure_20260928/ob_complete_strata.pdf),
[PNG](ob_complete_strata_figure_20260928/ob_complete_strata.png) and
[SVG](ob_complete_strata_figure_20260928/ob_complete_strata.svg) display all
eight retained methods in all fourteen frozen strata. Panels show weighted
F1, precision and recall on a common 0-100 percent scale, with fixed method
order. All 288 finite values and 48 missing cells are retained. Gray NA cells
are empty bins, not zero scores. Counts refer to reference families; strata
overlap across feature dimensions.

The source is the [independently checked descriptive export](OB_COMPLETE_STRATA_RESULT_20260928.md),
pinned by SHA-256. This figure adds no confidence intervals, significance
tests, causal interpretations or independent confirmation. It supplements,
not replaces, the original three-method uncertainty figure. Descriptors do
not establish fragmentation, domain structure or evolutionary duplication.
Full-reference scoring and low-certainty conventions remain unchanged.

The PNG was visually inspected: cell values, method labels, family counts,
color bar and notes fit without clipping or overlap. Automated rendering
checks compare every displayed cell and image-array value with the source,
including NA cells, and check text bounds. Output identities and the plotting
source are retained in the [manifest](ob_complete_strata_figure_20260928/manifest.json).

```sh
python -B -m benchmark_tools.plot_ob_complete_strata --source benchmark_tools/results/ob_complete_strata_20260928/report.json --output /tmp/ob-complete-strata-figure
```

The output directory must be unused. Reproduction uses the retained score
table, not raw inference. The main manuscript and portable figure bundle
have not yet been refreshed to include this supplementary figure.
