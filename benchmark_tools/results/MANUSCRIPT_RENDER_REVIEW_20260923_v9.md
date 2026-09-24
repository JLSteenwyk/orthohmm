# Working Manuscript Review, Version 9

[HTML](PUBLICATION_MANUSCRIPT_REVIEW_20260923_v9.html) and
[PDF](PUBLICATION_MANUSCRIPT_REVIEW_20260923_v9.pdf) keep the three supplementary
figures with their complete captions. Duplication and sequence-identity panels
and captions occupy page 22; the fragment panel and caption occupy page 23.
Both pages were visually checked at 1300-pixel maximum page dimension.

The manuscript differs from v8 only by three presentation-only fenced Divs.
Removing those wrappers produces a Pandoc AST exactly equal to the manuscript
at `aee085e`; no text, image, endpoint, scientific result or citation changed.
The [asset audit](manuscript_asset_review_20260923_v9.json) retains 198 local
occurrences and 180 targets. Fourteen focused renderer/artifact tests pass,
including a new HTML check for all three figure-caption groups.

The PDF remains 37 pages with 11 embedded figures and no detected out-of-page
text/image blocks. [Review receipt](manuscript_pdf_review_20260923_v9.json)
records output hashes and scope. Other pages were not re-reviewed visually for
v9; the full-page overview is retained in [v8](MANUSCRIPT_RENDER_REVIEW_20260923_v8.md).
Small figure labels, condensation, copyediting and journal-specific layout still
need attention. This remains an incomplete working publication package.
