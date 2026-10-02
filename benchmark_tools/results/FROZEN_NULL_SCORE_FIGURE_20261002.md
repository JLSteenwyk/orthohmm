# Synthetic Null Score Figure And Manuscript Integration

The [nine-panel figure](figures_frozen_null_scores_20261002_v2/frozen_null_scores.pdf)
and [all 90 plotted values](figures_frozen_null_scores_20261002_v2/plotted_values.tsv)
are generated from the retained synthetic observations. The current main and
extended manuscript sources now describe the protocol, outcome, bounded
interpretation and availability. This is not a calibration fix, new native
experiment or publication readiness. All biological scores/defaults remain frozen.

## Data And Display

The exporter independently verifies the supplied compressed-result and receipt
hashes, checks the decoded result against the original audit pin, and recomputes
all 90 endpoint counts/intervals from the saved raw scores. It rejects missing
or duplicated cells, changed statistics, negative scores, paired-band changes
and unsupported scope flags. This repeats numerical readback only, not native
scoring. Source and direct input/output pins are retained in the
[export manifest](figures_frozen_null_scores_20261002_v2/manifest.json).

All composition/length cells and cutoffs remain visible. Full/width-64 glyphs
are slightly horizontally offset; the exact cutoff coordinates remain in the
TSV. Whiskers are exact intervals with the full 90-endpoint adjustment.
Zero-hit triangles show upper limits, not positive observed frequencies. The
dashed model-tail reference is a diagnostic hypothesis, not a verified null
distribution. Forced-score tails are not real-data orthology error rates.

The first export required a presentation revision: length annotations were
too close to near-one tail glyphs in glutamine panels. Preserve that export
locally; move only the labels into clear lower-right space and generate a new
v2 directory. First and final plotted TSVs are byte-identical. The
[manual review receipt](frozen_null_score_figure_review_20261002.json) preserves
the earlier artifact hashes and the reason for revision rather than overwriting
the failed visual attempt.

Both the final PNG and a bitmap rendered from its one-page PDF were inspected.
All nine panels, labels, glyph meanings and notes are readable, with no observed
clipping or overlap. PDF block bounds have zero violations. The combined
plot/raw-readback/scoring/main-prose/render/print/PDF/artifact panel passes
**93 cases in 8.00s**, zero errors, failures or skips. The earlier 50-case
figure/scoring panel passed before the visual annotation issue was found;
passing tests did not replace visual inspection. A new overlap check now
ensures length labels avoid the actual marker bounding boxes.

## Reproduction

From the repository, use a new destination for every export:

```bash
python -m benchmark_tools.plot_frozen_null_scores \
  --scores benchmark_tools/results/frozen_null_score_observations_20261002.json.gz \
  --scores-sha 590db81716bc0af78acb66558fe0d0bb91356f99be90d229f1aa9bcdc3dc8c16 \
  --receipt benchmark_tools/results/frozen_null_score_receipt_20261002.json \
  --receipt-sha 7e059f8486f4801be7b2e5ed670464d0ac936111f45a76fe58a4458ec5b03ab9 \
  --output /absolute/new/null-score-figure
```

This does not load the native scorer or run inference. Output file bytes may
differ between rendering-library versions; displayed values and the explicit
inputs must remain identical. The one-page PDF is 30,735 bytes, SHA-256
`1ddd5fde7fd2855ad7f79e7f230a1c660e9b0fbe0a431e6bbb270eddf8131b01`.

## Remaining Review

The manuscript sources include this finding, but the earlier seven-page
HTML/PDF and review archives remain exact historical snapshots. A separately
dated main-text render and visual review are still pending at this source
checkpoint. No controlled timing, data-rights clearance, other-QfO uncertainty,
complete study restoration or public release follows from this figure.
