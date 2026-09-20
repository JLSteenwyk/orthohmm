# Explicit Manuscript Citation Inventory

The [initial inventory](manuscript_citation_inventory_20260920.json) of the
manuscript at `25b9615` matched 24 of the 37 selected v5 bibliography entries.
It found no unmatched DOI links, but 13 entries lacked an explicit matched
link. One of those, igraph, was cited through a documentation URL differing
from its bibliography URL; this was a matching limitation, not absence of
all attribution. Other references were principally present in supplements.

The Methods now places companion-tool, numerical-library and annotation-resource
references beside their documented roles. It uses the selected igraph URL
and retains the distinctions between method citations, invoked versions,
historical annotation snapshots and later resource descriptions. In particular,
the 2026 GO paper is not described as the source of QfO2020 annotations.

The [updated inventory](manuscript_citation_inventory_20260920_v2.json)
matches all 37 selected entries to 37 explicit citations, with no unresolved
DOI/citation-ID records or other unmatched external links. All original
manuscript evidence links remain separate local-link observations, not
bibliographic references or automatically validated evidence.

## Reproduction

```sh
python -m benchmark_tools.audit_manuscript_citations \
  --manuscript benchmark_tools/results/PUBLICATION_MANUSCRIPT_DRAFT_20260916.md \
  --bibliography benchmark_tools/results/publication_bibliography_20260920_v5.csl.json \
  --output NEW_CITATION_INVENTORY.json
python -m pytest -q tests/unit/test_audit_manuscript_citations.py
```

The auditor uses Pandoc's structured Markdown AST, not textual matching of
link-like text in code blocks. It inventories DOI resolver links (including
case/URL-encoding variants), exact bibliography URLs, and Pandoc citation IDs.
Duplicate IDs/DOIs are rejected; ambiguous shared URLs and unknown citation
IDs remain unresolved. Each occurrence carries its section and source record.
Report and source hashes establish which manuscript was checked. Later draft
changes require a fresh inventory; the earlier report is not silently updated.

Nine auditor tests and 25 bibliography/rendering tests pass (34 total).
Tests include actual Pandoc parsing, code blocks, reference-style links,
unknown identifiers, ambiguous URLs and DOI normalization. The data-derived
37/37 result is not a claim that every scientific statement is supported,
all used software is cited, sources were interpreted correctly, or bibliography
semantics/rights/journal formatting are complete. Those reviews remain open.
