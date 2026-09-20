# Working Manuscript Render Review

The [HTML working draft](PUBLICATION_MANUSCRIPT_REVIEW_20260920.html) renders
the manuscript with its repository-relative evidence links and eight inline
figures. Each figure is linked to its original full-resolution PNG, preserving
the alt text and surrounding scientific caption. This is a review artifact,
not submission-ready typesetting or a standalone portable archive.

An initial Chrome print of the pre-link-wrapper draft (`56ab41c`) produced
32 letter-size pages. PDF image inventory found all eight images. Pages
5, 11, 12, 14, 15, 18 and 21, which contain all inline figures, were visually
inspected: figures are visible without obvious clipping, but default-width
labels are small. The HTML links added afterward provide full-resolution
access without changing the scientific figure files. The final HTML link
structure is checked by HTML parsing; no full browser interaction or mobile
layout test is claimed. The other 25 initial PDF pages were not visually
reviewed, and the final linked HTML has not received a new complete print review.

The [local-asset inventory](manuscript_local_assets_review_20260920.json)
records 180 local link/image occurrences, 162 unique targets and their hashes.
All targets exist and were Git-tracked at review time; no missing target was
silently omitted. This does not validate fragment anchors, linked scientific
contents, transitive provenance or redistribution rights. All 37 selected
bibliography entries remain matched in the
[refreshed citation inventory](manuscript_citation_inventory_20260920_v3.json).

## Reproduction

Run from the repository root with a fresh sibling output name in
`benchmark_tools/results`, so relative links resolve correctly:

```sh
pandoc --from=markdown --to=html5 --standalone \
  --metadata=title:"OrthoHMM publication working draft" \
  benchmark_tools/results/PUBLICATION_MANUSCRIPT_DRAFT_20260916.md \
  --output benchmark_tools/results/NEW_MANUSCRIPT_REVIEW.html
google-chrome --headless --disable-gpu --no-pdf-header-footer \
  --print-to-pdf=ABSOLUTE_NEW_PDF \
  file:///ABSOLUTE_REPOSITORY/benchmark_tools/results/NEW_MANUSCRIPT_REVIEW.html
python -m pytest -q tests/unit/test_manuscript_review_artifacts.py
```

The reviewed print remains local at
`benchmarks/work/publication_manuscript_visual_20260920/manuscript.pdf`;
its hash and visual scope are in the inventory. PDF timestamps are not
byte-reproducible. This dated HTML remains fixed when the living draft changes;
future renders and source inventories must be separately identified. Final
figure sizing, manuscript condensation/copyediting, citation adequacy,
remaining results and journal-specific formatting remain open.
