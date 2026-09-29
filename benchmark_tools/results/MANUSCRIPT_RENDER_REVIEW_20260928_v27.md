# Extended Manuscript Review v27

The [HTML](publication_manuscript_review_20260928_v27.html) and
[47-page PDF](publication_manuscript_review_20260928_v27.pdf) now include
the complete eight-method OrthoBench descriptive strata figure and caption.
This also refreshes intervening source-text changes since the historical
v26 rendering; it is not a new scientific admission.

The [HTML asset receipt](manuscript_asset_review_20260928_v27.json) records
273 local references, 259 distinct tracked targets, no untracked targets
and empty Pandoc parse/render stderr. Chrome printed the local HTML with
a separate `/tmp/orthohmm-review-v27-profile` profile and no PDF headers or
footers. The process exited zero. After printing, all 265 HTML/source/target
records were rehashed alongside the PDF, asset receipt and layout checker.

The [layout receipt](manuscript_layout_review_20260928_v27.json) finds no
text/image block outside page bounds at one-point tolerance. The reusable
checker rasterized pages 18-21 around the required new paragraphs. Pages
19 and 20 were visually inspected: the earlier uncertainty panel, new
descriptive figure, captions and adjoining text have no observed clipping
or overlap. The all-method heatmap is dense at manuscript width; individual
cells should be consulted in the linked full-size figure or score table.
Pages 18 and 21 were not visually reviewed in this pass, nor were all other
pages. `visual_review_complete` therefore remains false. Sixteen renderer
and layout-checker tests pass, including missing-text and changed-source
rejection and adjacent-page selection.

```sh
python -B -m benchmark_tools.review_manuscript_pdf --pdf benchmark_tools/results/publication_manuscript_review_20260928_v27.pdf --assets benchmark_tools/results/manuscript_asset_review_20260928_v27.json --output /tmp/orthohmm-v27-layout-new --phrase 'The descriptive extension retains' --phrase 'Supplementary figure: all fourteen'
```

The output directory must be unused. Historical local-asset receipts bind
specific file versions; future source or linked-report edits may correctly
prevent their revalidation. No image rendering, text search or bounds check
establishes scientific correctness, full readability, rights clearance,
controlled timing or publication readiness. The concise main-text draft
and its separate six-page preview are unchanged.
