# Corrected-Score And Resource Manuscript Refresh

Updated `PUBLICATION_MANUSCRIPT_DRAFT_20260916.md` to include the first
independently admitted corrected reconciliation score, its six-endpoint
comparison with p0_c0_r0, precision-recall trade-offs, and the five-of-eight
factorial completion boundary. The draft no longer states that all R-on
assessments are unfinished. Neither p0_c0 arm is labeled satellite_v2.

Added the completed 18-task overhead audit and retained host-pressure
diagnostic to the resource limitations section. Measurement failures and
missing scheduler records remain separate from inference success; only one
paired comparison is available. The original 27 runs remain descriptive.
No timing, accuracy-significance, or publication-readiness claim was upgraded.

## Verification

- Pandoc 3.1.3 parsed the full GFM draft into JSON and rendered plain text.
- Parsed Link/Image nodes contain 129 local occurrences and 124 distinct
  targets; every target resolves to an existing file relative to the draft.
- Extracted all 12 numerical cells in the new reconciliation table from
  the Pandoc AST and compared them with the two source-manifest rows,
  rounded to six decimal places. All match.
- `git diff --check` passes for the manuscript changes. No scientific code
  changed, so native inference and unit tests were not rerun for this edit.

Verified manuscript SHA-256:
`c08ad1cfd8c083383a33c362faec523f8ee98567b0da5f182f859704ceb3ddbe`.
Source factorial manifest SHA-256:
`dcdca371ea006987383a7a06cd2178afb619d3c65fe73f6dd4eeb79a2eda07d6`.

The parse and rendering remain at
`benchmarks/work/publication_manuscript_ast_20260918.json` and
`benchmarks/work/publication_manuscript_plain_20260918.txt`. Regenerate the
parse from repository root with:

```bash
pandoc --from=gfm --to=json \
  benchmark_tools/results/PUBLICATION_MANUSCRIPT_DRAFT_20260916.md \
  --output=benchmarks/work/publication_manuscript_ast_20260918.json
```

This is a local-link and selected-table consistency check, not a verification
of every manuscript assertion, external URL, reference license, statistical
analysis or figure layout. The draft remains incomplete pending corrected
comparators/factorial, controlled resource evidence, and the other documented
scientific and release requirements.
