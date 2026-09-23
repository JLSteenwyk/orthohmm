# Refreshed Working Manuscript Review

The [23 September HTML](PUBLICATION_MANUSCRIPT_REVIEW_20260923.html)
includes the completed corrected FastOMA results, seven-method SwissTrees
uncertainty, descriptive secondary strata, and the explicit DGX deferral.
It links the newly bundled comparison PNG, byte-identical to the earlier
22098 PNG. The 20 September review is preserved unchanged.

The new render command uses Pandoc's parsed document to check local Link
and Image targets before writing HTML, refuses missing/escaping assets and
existing outputs, and records source, tool, target and output hashes. It
requires the HTML beside the manuscript so relative paths remain valid.
The [asset receipt](manuscript_local_assets_review_20260923.json) records
183 occurrences, 165 unique files, no untracked targets and no Pandoc
warnings. All eight HTML figures link to existing full-resolution images.
The [citation inventory](manuscript_citation_inventory_20260923.json)
matches all 37 selected bibliography entries, with no unresolved explicit
citations. This does not establish scientific citation adequacy.

Reproduce from the repository root with fresh sibling filenames:

```sh
python -m benchmark_tools.render_manuscript_review --repo . \
  --manuscript benchmark_tools/results/PUBLICATION_MANUSCRIPT_DRAFT_20260916.md \
  --output benchmark_tools/results/NEW_MANUSCRIPT_REVIEW.html \
  --report benchmark_tools/results/NEW_MANUSCRIPT_ASSETS.json
python -m benchmark_tools.audit_manuscript_citations \
  --manuscript benchmark_tools/results/PUBLICATION_MANUSCRIPT_DRAFT_20260916.md \
  --bibliography benchmark_tools/results/publication_bibliography_20260920_v5.csl.json \
  --output benchmark_tools/results/NEW_MANUSCRIPT_CITATIONS.json
python -m pytest -q tests/unit/test_render_manuscript_review.py \
  tests/unit/test_manuscript_review_artifacts.py \
  tests/unit/test_audit_manuscript_citations.py
```

The render and its links were structurally checked, not browser-reviewed or
printed to PDF this turn. External URLs, fragment anchors and embedded raw
HTML are not validated by the renderer. Target identities capture this
review's state; living documents may later change. This remains a dated
review requiring repository assets, not a standalone archive, complete
reproduction package, final typesetting or publication-readiness finding.
OrthoMCL recovery, queued robustness results and controlled resource evidence
remain incomplete. No DGX access occurred after the user's deferral.
