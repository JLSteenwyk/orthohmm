# Condensed Main Text Review

The [three-page PDF](publication_main_review_20260927_v2.pdf) and
[HTML](publication_main_review_20260927_v2.html) render the unchanged
[condensed manuscript](PUBLICATION_MAIN_TEXT_20260927.md). These are working
review documents, not a complete publication bundle or journal submission.

The initial four-page local render exposed a duplicated title and minimal
print-page margins. The renderer now uses Pandoc's page-title metadata without
inserting another visible H1. A hash-recorded print header specifies 18 mm
page margins, keeps headings with following text, and requests three-line
widow/orphan control. Historical renders were not overwritten. All eight
renderer tests pass, including a real Pandoc regression check for one visible
H1, retained HTML title, print rules and header provenance.

The [render receipt](publication_main_render_20260927_v2.json) records the
source, renderer, print header and 12 tracked local evidence targets. Source
and target hashes were rechecked after Chrome printed the PDF. The ledger
was subsequently updated to record this review; its receipt hash describes
the pre-update snapshot, not the evolving ledger forever. Links require the
repository layout; this is not a portable archive of all evidence.

The [layout receipt](publication_main_layout_20260927_v2.json) records three
pages and zero text/image-block bounds violations at one-point tolerance.
All three pages were subsequently rasterized and visually inspected at 1.2x:
text and headings are legible, margins are present, and no clipping or overlap
was observed. Paragraphs continuing across pages remain readable. The receipt's
`visual_review_complete=false` precedes inspection; this note records the
completed visual check of this three-page draft only. It does not claim review
of the extended manuscript, linked figures or a journal's layout requirements.

Reproduction uses `benchmark_tools.render_manuscript_review` with the condensed
Markdown, fresh sibling HTML and fresh report paths, followed by headless
Google Chrome with a separate temporary profile, `--no-pdf-header-footer`
and `--print-to-pdf`. No inference, scoring, endpoint, numerical claim or
scientific default changed. Uncertainty, TreeFam-source, controlled timing,
rights and archival release requirements remain open.
