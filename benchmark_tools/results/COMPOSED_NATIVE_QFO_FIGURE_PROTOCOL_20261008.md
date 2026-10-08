# Composed-Native QfO Figure

## Prerequisites

This is a prospective figure/table consumer, not a final benchmark figure.
Require the actual composed five-admitted-cell snapshot and native-count paired
bootstrap, with successful final native12review/conversion/six-endpoint scoring/
independent admission and raw-family audit. Do not render from synthetic test
data or the five completed but unadmitted endpoints of failed24038scoring.
If final12cannot be admitted, preserve the existing four-cell figure and report
the actual missing outcome in a separate addendum rather than manufacturing
a fifth cell or repeating inference.

Preserve the frozen four-cell plot, reviewer, assets and manuscript. The new
sources are `plot_composed_native_qfo_scores.py` and
`review_composed_native_qfo_figure.py`.

## Scientific Content

The plot consumes the explicit composed snapshot and newly computed native
intervals, not an ordinary-schema translation or cached uncertainty substitution.
It replays the scientific reporting snapshot in the retained scientific Python
environment, checks input/source/evidence hashes before and after rendering,
and independently checks F1 against native precision/recall, coverage counts,
interval definitions, arithmetic and missing-cell scope.

- All seven cell statuses remain in `status.tsv`.
- `scores.tsv` has42rows:30admitted values and12missing values, not zeros.
- `coverage.tsv` has seven rows, including genuine failed-scoring coverage
  without accuracy admission. Missing inference coverage remains missing.
- PanelsA/B plot15orthology F1 scores and15GO/EC/FAS similarities; GO, EC and
  FAS are not F1. The secondary six-metric mean is not plotted.
- PanelC plots the five admitted cells' all-input relation coverage, not
  reference-gene accuracy. Failed-scoring coverage remains in the table.
- PanelD plots only estimable conditional SwissTrees differences, with nominal
  and42endpoint-adjusted intervals from100000new shared-family draws. It does
  not assume sign, zero inclusion or superiority. All14contrast definitions
  and missing required cells remain in the input bootstrap report.
- P denotes profile refinement, C candidate expansion and R phylogenetic pair
  inference; initial HMM search remains on in every cell. R also changes
  clique-derived versus inferred-pair prediction semantics.

Recovered failed timing remains ineligible. No timing is plotted or isolated
efficiency claimed; the complete shared-Threadripper contention disclosure is
retained in the figure manifest. The18development-exposed SwissTrees families
retain exchangeability, selection and approximate percentile-coverage limits.
This is not independent confirmation, raw admission or a full-factorial claim.

## Commands After Successful Admission

Use the figure Python environment with Matplotlib, NumPy, Pillow and PyMuPDF;
set `SCIENTIFIC_PYTHON` to the retained Python3.10scientific environment used
for snapshot replay. Disable user-site/Python/dynamic-loader injection and use
one BLAS/OpenMP thread. Use actual artifact paths/digests and fresh output paths.

```bash
python -B -m benchmark_tools.plot_composed_native_qfo_scores \
  --snapshot "$SNAPSHOT" --snapshot-sha256 "$SNAPSHOT_SHA" \
  --native-intervals "$NATIVE_INTERVALS" --native-intervals-sha256 "$NATIVE_INTERVALS_SHA" \
  --validation-python "$SCIENTIFIC_PYTHON" --output "$FRESH_FIGURE_ROOT"

python -B -m benchmark_tools.review_composed_native_qfo_figure \
  --figure-root "$FRESH_FIGURE_ROOT" --pdf-preview "$FRESH_PREVIEW" \
  --output "$FRESH_READBACK"
```

The reviewer independently reads all table values against the bound snapshot
and intervals, verifies null entries and statuses, decodes PNG/PDF, checks all
five cell colors, required SVG/PDF labels, page count/dimensions and PDF text
bounds. It does not call the plot's data extractor to establish table equality.
Machine checks leave `visual_review_complete=false`; manually inspect actual
production assets before publication integration. Figure generation makes zero
bootstrap draws but explicitly consumes the prior100000new native-count draws.

## Test Evidence

Test handoffs explicitly invent a final cell and stub scientific replay. Real
count/statistic/bootstrap, figure, CSV and PDF/PNG/SVG kernels run unchanged.
Tests reject partial admissions, failed-score imputation, wrong endpoint types,
F1 arithmetic, coverage denominators, source/cross-snapshot bindings, changed
bootstrap scope/contrast definitions, invalid interval nesting and resealed
table/asset mutations. Positive intervals need not include zero. Test assets
are not actual native12outputs and must not be incorporated as scientific figures.
