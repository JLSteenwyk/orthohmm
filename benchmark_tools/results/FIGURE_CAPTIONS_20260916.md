# Audited Accuracy Figures

These figures describe development-exposed historical results, not a frozen
publication baseline or an independent superiority test. Source values are
`publication_comparison_20260916.json`; each plotted coordinate and interval,
the plotting source hash, input hash, software version, and artifact hashes
are recorded in `figures_accuracy_20260916/manifest.json`.

## Accuracy Overview

`accuracy_overview`: (A) OrthoBench weighted precision and recall for the
eight retained method configurations. (B) Three Kingdoms micro pair F1 in
the BUSCO-reference gene universe. The two panels evaluate different
relationships and must not be averaged or treated as replicates. BUSCO-only
scoring excludes false positives involving genes outside the reference
universe. Stable symbols/colors identify methods; no rank-based selection
or truncated score axis is used. OrthoFinder sequence denotes its MCL
checkpoint, and FastOMA used a supplied tree. Competitor raw-output
provenance gaps described in the comparison report remain open.

## QfO Endpoints

`qfo_endpoints`: native precision/recall coordinates for VGNC, SwissTrees,
and TreeFam-A, and functional similarity versus assessed relation count for
EC, GO, and FAS. Relation counts are not total prediction coverage. The
functional panels use logarithmic count axes; all similarity/precision/
recall axes span zero to one. No custom six-metric mean is plotted.
Recorded uncertainty fields have heterogeneous definitions across native
implementations and are not presented as a common standard error or CI.
OrthoMCL final-group scoring remains pending in this source snapshot and is
not plotted; its pre-MCL graph is not substituted. The legend retains its
identity so the same method encoding can be used in the completed update.

## Paired Differences

`orthobench_paired_differences`: OrthoHMM minus OrthoFinder 3.1.5 full in
weighted F1, precision, and recall, in percentage points. Points are observed
differences; thick segments are nominal 95% paired percentile intervals,
and thin segments use Bonferroni-adjusted tails across six contrasts/metrics.
The statistic is recomputed from sufficient counts for 20,000 paired
reference-family bootstrap replicates (seed 20260916), not averaged from
per-family F1. Intervals assume approximately exchangeable reference
families and do not account for prior method selection on OrthoBench.
Both F1 intervals include zero; the plot does not establish an F1 advantage.

## Reproduction

From the repository root, use a fresh output directory:

```bash
python benchmark_tools/plot_publication_accuracy.py \
  --comparison benchmark_tools/results/publication_comparison_20260916.json \
  --output benchmark_tools/results/figures_accuracy_reproduced
```

PNG (200 dpi), PDF, and SVG exports are generated together. The script
refuses existing output directories, rejects inconsistent axes and paired
estimates, and omits all QfO coordinates unless the method's recorded status
is `metrics_available`. Regenerate the comparison after verified OrthoMCL
workflow completion, then generate a new versioned figure bundle rather
than overwriting this pending-results snapshot.

This is an accuracy-figure milestone only. Independent validation, the
method diagram, ablation figures, error strata, scaling, biological
application, manuscript text, and release package remain unfinished.
