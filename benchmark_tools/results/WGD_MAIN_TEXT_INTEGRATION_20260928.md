# Quantitative Biological Application In Main Text

The main manuscript now reports the prespecified duplicate-gene application's
counts, denominators, homolog coverage and adjusted paired contrasts rather
than only describing the comparison qualitatively. No method was changed,
new biological prediction generated, example replaced or endpoint retuned.

Among 239 input-eligible pairs, phylogenetic OrthoHMM separated 238, high
sensitivity 58 and full OrthoFinder 236. Among 231 reference-eligible pairs,
supported separations were 193, 56 and 227, respectively. Mean homolog coverage
was 82.338%, 99.149% and 98.413%. Coverage uses the union of anchor groups:
merged anchors can retain high coverage, and this is not orthology recall.
The text retains the negative comparisons with full OrthoFinder and links
all five methods rather than selecting only favorable results.

The existing 20,000-resample, 12-endpoint-adjusted contrasts are reproduced,
not newly fitted. Each eligible pair occupies a distinct reference pillar;
conditional pillar uncertainty is not evidence of independent generalization
or cross-species copy-specific orthology. Six prespecified case traces remain
linked. Five focal homologs leave anchor groups during root-lineage
reconstruction before satellite constraints, but remain elsewhere in the
output. The trace does not distinguish topology from lineage-rule causation.

## Verification

- [Native-output audit](biological_wgd_main_text_recheck_20260928.json)
  reproduced all five methods, all 240 pair records, six examples and twelve
  contrasts. A resume-time rerun reproduced this receipt exactly.
- Regenerated case traces are byte-identical to the retained historical trace;
  the duplicate stays local rather than entering the repository again.
- All 51 focused audit, case-trace, scoring, bootstrap and renderer tests pass.
- [Render receipt](publication_main_render_20260928_v7.json) records 16
  citations and 16 local link occurrences across 15 tracked targets, with no
  Pandoc warnings. All five PDF pages were visually inspected without clipping
  or incoherent overlap; numerical labels and word bounds were checked.
- [Layout and evidence receipt](publication_main_layout_20260928_v7.json)
  binds the corrected manuscript, source results, native audit and PDF.

Read the [five-page review PDF](publication_main_review_20260928_v7.pdf) or
[HTML export](publication_main_review_20260928_v7.html). Earlier exports remain
historical; local v6 was an intermediate wording revision, not the current
manuscript. This is scientific reporting and reproduction, not new validation.
Controlled timing, remaining uncertainty work and release requirements are
still open. No benchmark was launched, unrelated job changed or DGX accessed.
