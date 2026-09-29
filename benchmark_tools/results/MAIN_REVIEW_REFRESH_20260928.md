# Main Manuscript Review Refresh

Generated [HTML](publication_main_review_20260928_v11.html) and
[PDF](publication_main_review_20260928_v11.pdf) from the current condensed
manuscript, preserving the scientific source and older review versions.
[Render receipt](publication_main_render_20260928_v11.json) binds source,
bibliography and direct local targets at render time.
[PDF check](publication_main_pdf_review_20260928_v11/report.json) records all
six rasterized pages and zero page-bound violations. All six images were
visually inspected: no clipping or incoherent overlap observed. The BUSCO
bibliography entry continues across pages 5-6; journal formatting remains open.
The automated receipt deliberately does not assert manual visual review.

Page 2 includes the complete 21-endpoint OrthoBench comparison: phylogenetic
OrthoHMM minus full OrthoFinder F1 +1.370 percentage points, adjusted interval
[-7.281, 12.166], precision +15.705 [1.175, 30.634], recall -13.151
[-26.947, -0.609]. These support a precision-recall trade-off, not established
F1 superiority. This rendering check does not independently validate statistics.

Seventeen focused renderer, PDF-review and complete-OrthoBench manuscript tests
pass. PDF generated with local headless Google Chrome, no print headers or
footers, from the rendered local HTML. No benchmarks were rerun. Slurm queue
was empty when inspected; this does not establish non-Slurm isolation.

Provenance exception: the initial v6 HTML render refused an existing path.
A subsequent print command nevertheless regenerated the pre-existing untracked
v6 PDF from its existing HTML. Original v6 PDF bytes were not preserved by that
command, so historical v6 PDF hashes must not be treated as current. No v6
files are included in this milestone. v11 was rendered successfully to fresh
paths, and its PDF print was guarded against overwriting an existing file.
Versions v7-v10 were left unchanged.

The linked progress ledger is mutable and changes after rendering; its receipt
hash is a snapshot, not a promise that the live ledger remains unchanged.
This is a review export, not a standalone archival bundle or publication-ready
release. Controlled timing, several QfO uncertainty estimates and the remaining
publication requirements remain open.
