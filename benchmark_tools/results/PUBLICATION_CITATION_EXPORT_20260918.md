# Selected Publication Citation Export

Exported 14 selected references from Crossref DOI metadata as CSL-JSON,
including complete deposited author lists. The explicit DOI/year selection
is `publication_citation_selection_20260918.json`; scientific source review
remains in `PUBLICATION_REFERENCES_20260917.md`.

The export retains the OrthoFinder 2026 correction separately and labels
the OrthoHMM 2024 record as a preprint, not a journal article. Crossref
`article-number` maps to CSL `number`, separately from actual `page` ranges.
The [CSL specification](https://docs.citationstyles.org/en/stable/specification.html#appendix-iv-variables)
defines number as an item identifier and page as a page range; it does not
define an article-number variable. Final journal-style rendering remains
to be checked, particularly styles that omit number for journal articles.

## Reproduction

Raw responses are retained under
`benchmarks/work/publication_citations_20260918/`, not committed.
The provenance report records source URLs, retrieval times, byte sizes and
SHA-256 values. This offline replay preserves original retrieval times,
rejects altered responses or selections, and requires new output paths:

```bash
python benchmark_tools/export_publication_citations.py \
  --manifest benchmark_tools/results/publication_citation_selection_20260918.json \
  --raw-directory benchmarks/work/publication_citations_20260918 \
  --cached-provenance benchmark_tools/results/publication_citation_provenance_20260918.json \
  --output /tmp/orthohmm-references.csl.json \
  --provenance /tmp/orthohmm-reference-replay.json
```

For a fresh acquisition, omit `--cached-provenance` and provide a new raw
directory. New metadata may differ; it is not automatically the frozen
snapshot. Raw responses must be included in the eventual archival bundle
after reuse-license review, or acquired and checksum-verified independently.

CSL output and independently replayed output both have SHA-256
`c44a7eb5b6eec9e4fa5d0c13146d194fe932888931fd8adb62fa66095f45f3e8`.
The committed provenance report has SHA-256
`4fc05a637ffac68cac37bbd5da7b39b9597b5977a01219534d101a8e5b4ecc6a`.
Sixteen focused tests pass, covering title markup, consortium names,
suffix preservation, preprint typing, identifiers, invalid metadata,
offline replay, checksum/inventory/URL mismatches and overwrite refusal.

## Remaining Limits

- Metadata is not full-text scientific validation, software-version
  provenance, or a license audit. Abstracts and cited-reference lists are
  excluded from the committed export.
- The OrthoHMM record deposits Thomas J. Buida without the III suffix;
  preserve the source metadata and resolve this explicitly before final
  submission rather than silently inventing a correction.
- This selection is not the complete dependency, simulator and dataset
  bibliography. Those citations and journal-specific rendering remain open.
- The TreeFam family-level source files remain unretrieved. A bibliographic
  record is not a substitute for missing reference-generation inputs.
