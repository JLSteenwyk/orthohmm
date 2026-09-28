# YGOB Overlap-Stratum Figure

[PDF](figures_ygob_overlap_20260928/ygob_overlap_strata.pdf),
[PNG](figures_ygob_overlap_20260928/ygob_overlap_strata.png),
[SVG](figures_ygob_overlap_20260928/ygob_overlap_strata.svg), and
[plotted values](figures_ygob_overlap_20260928/plotted_values.tsv).
The [manifest](figures_ygob_overlap_20260928/manifest.json) records source,
plotter, Matplotlib version and generated-file identities.

## Caption

Descriptive YGOB curated-group recovery partitioned by the retained
homology-overlap screen. Panels A-C show screen-positive pillars; D-F show
pillars without a qualifying hit. All four frozen methods and all three
metrics are retained on common 0-100% scales. Points are ratios of summed
original TP/FP/FN statistics, not averages of pillar scores. Original
cross-pillar false-positive allocations, including penalties across strata,
remain included. Predictions were not filtered and rescored. Counts above
each row show the different reference composition and singleton prevalence.
The OrthoFinder sequence-only checkpoint is diagnostic, not a separate
finalized phylogenetic analysis. These co-membership scores include
within-species pairs and do not represent resolved pairwise orthology.

No confidence intervals or significance claims are added. No-hit status does
not establish absence of remote homology or independent annotation ancestry.
The screen-negative F1 advantage of phylogenetic OrthoHMM accompanies lower
recall; it is not evidence of better remote-ortholog recovery or a causal
effect of development exposure. See the
[prespecified descriptive protocol](YGOB_OVERLAP_STRATA_PROTOCOL_20260928.md)
and [analysis result](YGOB_OVERLAP_STRATA_RESULT_20260928.md).

## Reproduction

```bash
python -B -m benchmark_tools.plot_ygob_overlap_strata \
  --results benchmark_tools/results/ygob_overlap_strata_20260928.json \
  --sha256 464d5d72a6d16c27b68d2a1204ee81bb118e99430b9b8a4d75c5c24942895b8c \
  --output /tmp/ygob-overlap-figure-fresh
```

The output directory must not exist. The plotter checks every plotted ratio
against exact-rational count arithmetic, requires the complete method/stratum
panel and consistent reference sizes, and rejects unsupported confirmation
or CI claims. Undefined ratios are labelled rather than drawn as measured
zeroes. Its checks validate rendering inputs, not raw prediction generation;
the original aggregation and scoring audits remain the upstream evidence.

The main-text draft now reports both strata and both directions of the
precision-recall trade-off. Historical exports remain retained and are not
silently overwritten. No source prediction, primary score, scientific
configuration or independence designation changed.

## Review

46 focused plotting, aggregation and group-scoring tests pass. One initial
test incorrectly counted TSV physical lines despite a quoted multiline label;
it now uses the CSV parser and verifies all 24 rows and all four labels.
The PNG and independently rasterized vector PDF were visually inspected.
All 24 numeric labels were found in PDF text, and extracted word bounds stay
within the page. This verifies presentation, not scientific independence.

The [refreshed five-page main PDF](publication_main_review_20260928_v5.pdf)
and [HTML](publication_main_review_20260928_v5.html) include the new paragraph
and figure link. All five pages were visually inspected; 16 citations resolve,
14 local-link occurrences have tracked targets, and Pandoc reported no warnings.
The [layout receipt](publication_main_layout_20260928_v5.json) retains checks
and limitations, including the bibliography entry continuing across pages
4-5 and the mutable progress-ledger link. Final journal formatting and the
publication-readiness requirements remain unfinished.
