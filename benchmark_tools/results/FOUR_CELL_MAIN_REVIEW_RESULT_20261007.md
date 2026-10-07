# Four Cell Main Manuscript Render And Review

The [four-cell main text](PUBLICATION_MAIN_TEXT_20261007_v4.md) now has a
[checked HTML review](PUBLICATION_MAIN_TEXT_20261007_v4_checked.html) and
[21-page PDF](publication_main_review_20261007_v4_checked_print/document.pdf).
All pages were visually inspected; numerical scope, missing outcomes and
negative results remain explicit. This is a bounded review of the exact new
source, not publication readiness or a complete executable study archive.

## Actual Execution

The [first bibliography-selection failure](four_cell_main_render_selection_failure_20261007_v1.json)
remains failed. Prospective correction `a0cf2805` was committed/pushed before
one render with the actual parent-bound October bibliography. No original
source, bibliography or renderer was edited. Exact successful commands below
used the same clean environment and GNU-time format recorded in the
[generation result](FOUR_CELL_MAIN_TEXT_RESULT_20261007.md), with the rendering
interpreter `benchmarks/work/release_alert_refresh_20261001/venv/bin/python`.

```bash
python -B -m benchmark_tools.render_manuscript_review \
  --repo . --manuscript benchmark_tools/results/PUBLICATION_MAIN_TEXT_20261007_v4.md \
  --output benchmark_tools/results/PUBLICATION_MAIN_TEXT_20261007_v4_checked.html \
  --report benchmark_tools/results/publication_main_review_20261007_v4_checked_assets.json \
  --bibliography benchmark_tools/results/publication_bibliography_20261007_v1.csl.json

python -B -m benchmark_tools.print_manuscript_review \
  --assets benchmark_tools/results/publication_main_review_20261007_v4_checked_assets.json \
  --output benchmark_tools/results/publication_main_review_20261007_v4_checked_print \
  --browser /opt/google/chrome/google-chrome

python -B -m benchmark_tools.review_manuscript_pdf \
  --pdf benchmark_tools/results/publication_main_review_20261007_v4_checked_print/document.pdf \
  --assets benchmark_tools/results/publication_main_review_20261007_v4_checked_assets.json \
  --output benchmark_tools/results/publication_main_review_20261007_v4_checked_pdf_review \
  --phrase '.' --phrase 'Four of seven' --phrase 'P1/C0/R1' \
  --phrase 'candidate separation' \
  --phrase 'Timing measurements were collected on a shared Threadripper'
```

All three actual exits were 0. GNU-time observations, not native tool timings:

| Step | Elapsed seconds | Maximum RSS KiB | Result |
| --- | ---: | ---: | --- |
| Corrected HTML render | 1.29 | 129,024 | Empty parse/render stderr; 19 citation IDs |
| Browser print | 1.70 | 196,972 | 21-page PDF, 391,777 bytes |
| PDF bounds and page render | 1.81 | 63,448 | Zero bounds violations; all 21 pages rendered |

The actual receipt chain was independently checked: source/assets/HTML match
print inputs, printed PDF matches the reader input, and all 21 decoded image
records still match their hashes. All 24 score strings also occur on PDF page
8. Actual text selectors locate the four-cell section on page 8, candidate
separation on page 9 and timing disclosure on page 5. The separate
[full visual review](FOUR_CELL_MAIN_VISUAL_REVIEW_20261007.md) records actual
observations rather than rewriting immutable automated receipts.

## Bound Outputs

| Artifact | SHA256 |
| --- | --- |
| HTML | `3501e87369e98324adf93ece959340f6b1d63f064c97de7f73c6ee5e83d4695d` |
| Asset receipt | `995f707cf3dcf43b73bc6f95a8fbdfb12decb10884e0ed11142a650e7927fb21` |
| Print receipt | `e4abfea87036c3e4252267603bcf7166d766ebfbb13f5c712e23d4dd67f0d681` |
| Printed PDF | `2af09ee94d1a11393de3ceb4b45c7a021e122c41421e40959df877d55fa9692d` |
| Bounds and page receipt | `88e054a69db01f6ae73e06ed78c641d6036d5d78a0bc7a136990d37ed9de5a8d` |

The [asset receipt](publication_main_review_20261007_v4_checked_assets.json)
checks 110 local link occurrences across 107 unique tracked targets and all
19 citation IDs. It is repository-relative review evidence, not a portable
transitive archive. Browser profile files are operational state, not review
deliverables. Old manuscripts, reviews, figures and archives retain their
original bytes and scope. No new inference, scoring, bootstrap, admission,
default selection, independent confirmation or DOI follows from this review.
