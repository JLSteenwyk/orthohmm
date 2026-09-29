# Main Manuscript PDF Refresh

The [six-page PDF](publication_main_print_20260929_v2/document.pdf) now
includes the manuscript's corrected VGNC deletion sensitivity paragraph and
complete OrthoBench interval discussion. This closes the PDF-review gap left
by the 29 September v1 HTML export. No manuscript prose or scientific
results were changed in this refresh.

The [fresh HTML receipt](publication_main_render_20260929_v2.json) checked
16 citation IDs and 23 local-link occurrences covering 22 tracked targets.
The checked printer verified direct input hashes before and after printing;
its [receipt](publication_main_print_20260929_v2/print.json) records a
117,690-byte PDF, SHA256
`e62be33f1ea841ff72747f421feae990a64a0d387ddc746926a5ef9ad71913e5`.
The browser sandbox was explicitly disabled for this trusted local input,
as in the earlier checked printing experiment. Browser profile/cache files
are not release artifacts.

The [bounds check](publication_main_pdf_review_20260929_v2/report.json)
found no out-of-page blocks at its one-point tolerance. Selectors for
Abstract, 16,844, 2,893 and References plus adjacent pages covered all six
pages. All six generated page images were visually inspected: no clipping,
overlapping text or blank pages were observed. The FastOMA reference begins
on page 5 and continues on page 6; this is readable draft pagination, not
final journal typography. Figures remain links rather than embedded panels.
The mechanical receipt correctly does not itself certify visual inspection.

## Executed Workflow

From the repository root, using existing Python and Pandoc installations:

```bash
python -B -m benchmark_tools.render_manuscript_review --repo . \
  --manuscript benchmark_tools/results/PUBLICATION_MAIN_TEXT_20260927.md \
  --bibliography benchmark_tools/results/publication_bibliography_20260920_v5.csl.json \
  --output benchmark_tools/results/publication_main_review_20260929_v2.html \
  --report benchmark_tools/results/publication_main_render_20260929_v2.json
python -B -m benchmark_tools.print_manuscript_review \
  --assets benchmark_tools/results/publication_main_render_20260929_v2.json \
  --output benchmark_tools/results/publication_main_print_20260929_v2 \
  --browser /usr/bin/google-chrome --no-sandbox
python -B -m benchmark_tools.review_manuscript_pdf \
  --pdf benchmark_tools/results/publication_main_print_20260929_v2/document.pdf \
  --assets benchmark_tools/results/publication_main_render_20260929_v2.json \
  --output benchmark_tools/results/publication_main_pdf_review_20260929_v2 \
  --phrase Abstract --phrase 16,844 --phrase References --phrase 2,893
```

Existing outputs must not be overwritten. Future executions need fresh paths
and a fresh render receipt: the linked progress ledger changes after this
review, so its retained hash is a historical snapshot. This is not a
standalone archive, transitive source verification, publication readiness
or validation of the manuscript's scientific claims.
