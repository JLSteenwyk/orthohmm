# Main Text Citation Integration

The [cited HTML draft](publication_main_review_20260927_v3.html) adds 16 formal
references from the existing reviewed 37-record CSL bibliography. It attributes
the OrthoHMM preprint lineage, clustering methods, retained comparators
(including the OrthoFinder correction), benchmark resources, DIAMOND and the
experimental whole-genome-duplicate source. It does not claim the papers prove
the versions, settings or outputs of our runs. Local execution-evidence links
remain separate. No new literature metadata was fetched or inferred.

The renderer's optional `--bibliography` argument pins the supplied CSL file
alongside the manuscript, renderer and print header. Parsed citation IDs must
be present and bibliography IDs must be unique nonempty strings. Citations
without an explicit bibliography are rejected, as are paths resolving outside
the repository. Existing no-citation manuscripts remain supported. Validation
failures precede writing HTML or report files.

Thirteen focused tests pass, including actual Pandoc citation rendering and
missing-ID, duplicate-ID, absent-bibliography and escaping-path rejection.
The [render receipt](publication_main_render_20260927_v3.json) has empty parser
and renderer stderr and no untracked local evidence targets. Independent HTML
parsing found exactly the 16 reference entries named by the manuscript's
parsed citation set. This verifies reference resolution, not scientific
adequacy of every attribution or complete dependency citation coverage.

Reproduce with fresh output and report paths:

```bash
python -m benchmark_tools.render_manuscript_review \
  --repo . \
  --manuscript benchmark_tools/results/PUBLICATION_MAIN_TEXT_20260927.md \
  --bibliography benchmark_tools/results/publication_bibliography_20260920_v5.csl.json \
  --output benchmark_tools/results/FRESH_MAIN_REVIEW.html \
  --report benchmark_tools/results/FRESH_MAIN_RENDER.json
```

Pandoc's default citation style is used; no journal-specific style is selected.
The previous three-page PDF predates these citations and remains unchanged.
The cited HTML has not received a new browser/PDF visual review. Numerical
results and scientific defaults are unchanged. The publication ledger was
updated after rendering; its pinned target hash represents that earlier
snapshot. Scientific uncertainty, controlled timing and archival release
requirements remain open.
