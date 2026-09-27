# Integrated Result PDF Review

The [v26 PDF](publication_manuscript_review_20260927_v26.pdf) renders the
unchanged [v26 HTML](publication_manuscript_review_20260927_v26.html), including
the independently admitted full OrthoBench result and reconstructed-base
fixture evidence. It complements the earlier
[HTML-only review](MANUSCRIPT_RENDER_REVIEW_20260927_v26.md).

Headless Google Chrome printed the local HTML with a separate profile and
`--no-pdf-header-footer`. All source, renderer, HTML and 248 local target
hashes were rechecked against the HTML asset receipt after printing.
The [layout receipt](manuscript_layout_review_20260927_v26.json) records
46 pages and no text/image block bounds violations at one-point tolerance.
The PDF is 3,424,844 bytes, SHA256
`bbff2a11072f5d8f27de5665baed78636bcc2de14942b049b23f644536509563`.

Pages 44 and 45 were visually inspected at 1.3x raster resolution. The new
admission and runtime-reconstruction paragraphs, including the latter's
page break, are legible without clipping or overlap. A text search also
selected pages 27 and 28 for rasterization; those pages were not visually
reviewed in this pass. `visual_review_complete` remains false because this
was not a full-document visual or journal-typesetting review.

Earlier snapshots remain unchanged. This review does not establish new
scientific findings, controlled timing, rights clearance, complete runtime
portability or publication readiness.
