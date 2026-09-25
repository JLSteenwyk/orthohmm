# Working Manuscript Review, Version 10

[HTML](PUBLICATION_MANUSCRIPT_REVIEW_20260925_v10.html) and
[PDF](PUBLICATION_MANUSCRIPT_REVIEW_20260925_v10.pdf) shorten the opening
by moving bibliography audit history into a
[linked provenance note](PUBLICATION_MANUSCRIPT_BIBLIOGRAPHY_PROVENANCE_20260925.md).
The original history beginning with "An initial" is preserved verbatim in
that note. Everything from `## Study Objective` onward matches commit
`f3d6295` byte-for-byte; the body after that heading has SHA256
`25fc1fc8eaecb4d7c3435bc63e0454effaa383d33ce628857a3998bf7b926746`.
No scientific text, figure, result or in-body citation changed.

The [asset receipt](manuscript_asset_review_20260925_v10.json) verifies
189 local occurrences and 175 targets, with no untracked targets. The
separate provenance note has 11 local links to 10 existing targets, checked
with the same Pandoc/local-assets parser. These are existence checks, not
transitive provenance, scientific-content or redistribution-rights audits.

Rendered with `google-chrome --headless --disable-gpu
--no-pdf-header-footer --print-to-pdf=OUTPUT file:///ABSOLUTE_HTML_PATH`.
PDF SHA256: `ee70c37471bc07d5b78fc022d229e4e6e7e9741e89f65442b61b646e71e39bbc`.
HTML SHA256: `47089228ff750868697919b6eb13388d24a260fcb0ffbda166afb77f3143c294`.
PyMuPDF inspection finds 37 pages, 11 image occurrences and no text/image
blocks outside page bounds with one-point tolerance.

Visually inspected page 1 at 1600-pixel maximum dimension and pages 21-22
at 1300 pixels. The opening now reaches Methods; all three supplementary
figures retain their complete captions on the same page. Other pages were
not visually re-reviewed after pagination changed. Small labels, broader
condensation, copyediting and journal formatting remain unfinished. Earlier
reviews are retained. This is not a publication-readiness declaration.
