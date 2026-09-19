# Consolidated Selected Bibliography

## Current 37-Record Export

`publication_bibliography_20260919_v2.csl.json` adds the reviewed igraph
article citation to the original 36-record selection. All original records
and the eight earlier source selections are unchanged. A ninth source,
`publication_igraph_reference_20260919.csl.json`, is an explicit manual
transcription, not a fabricated Crossref export. The v1 files remain intact.

The [official Python citation guidance](https://python.igraph.org/en/0.11.6/)
recommends the 2006 article. The [official author/citation page](https://r.igraph.org/authors.html)
provides its BibTeX fields and author diacritics. Only the 2006 article was
selected, not the adjacent R-package citation or 2023 preprint. The official
BibTeX `pages = 1695` is preserved as CSL `page`; this is source transcription,
not independent verification of pagination. No article DOI is asserted.

`publication_igraph_reference_provenance_20260919.json` records these decisions
and checksums of both locally retained HTML snapshots. The new snapshot is
`benchmarks/work/igraph_authors_guidance_20260919.html`, 16,586 bytes,
SHA-256 `d37a52c9529431105d15f5c508969ef0c9114df9fc1c2280d4f69845cd6ae68a`.
Raw webpages are not committed.

```sh
python benchmark_tools/assemble_publication_citations.py \
  --manifest benchmark_tools/results/publication_bibliography_selection_20260919_v2.json \
  --output NEW_BIBLIOGRAPHY.csl.json --provenance NEW_PROVENANCE.json
```

The new output SHA-256 is
`a6f1c52694cd892f76dab1e940438540029b803b27b8296f0448fa714f0fd6f4`.
Thirty-two focused tests pass, including relocation of all nine actual source
exports, unchanged original records, exact igraph author/year fields, and
absence of a substituted DOI. The original output still has SHA-256
`6ca17ec5bcf85b24386db04f3c2abcf62d15cc296b130c5319905392f452a87f`.

This resolves only the igraph CSL omission. The OrthoHMM author suffix,
journal-specific rendering, full citation coverage and rights review remain
open. No scientific executor, benchmark result or default parameter changed.

## Original 36-Record Export

`publication_bibliography_20260919.csl.json` combines 36 selected records from
eight explicitly selected exports. It is a single bibliography input, not a
claim that every manuscript dependency has been cited or journal-formatted.
The assembly manifest pins every source filename and SHA-256; the provenance
report maps each citation ID to its input export and identifies the assembler.

The reviewed resource and QfO-service byline exports are selected instead of
their raw Crossref counterparts. The TreeFam paper retains 15 authors and
the 2022 QfO paper retains 31 top-level byline entries. No source fields are
rewritten: article identifiers remain distinct from pagination, website
access dates are not issue dates, the OrthoHMM record remains a preprint,
and the OrthoFinder 2026 correction remains a separate citation.

## Reproduction

Run from the repository root, using new output paths:

```sh
python benchmark_tools/assemble_publication_citations.py \
  --manifest benchmark_tools/results/publication_bibliography_selection_20260919.json \
  --output NEW_BIBLIOGRAPHY.csl.json --provenance NEW_PROVENANCE.json
```

No network request is made. The assembler rejects changed source bytes,
repeated source files, duplicate citation IDs or case-insensitive duplicate
DOIs. It preserves input and record order; it does not silently deduplicate
or choose among conflicting metadata. The source exports' original acquisition
and byline-correction provenance remain separate retained files.

The retained output and independent second assembly both have SHA-256
`6ca17ec5bcf85b24386db04f3c2abcf62d15cc296b130c5319905392f452a87f`.
All 31 focused assembly/export tests pass, including a relocated copy of the
actual selection and eight sources, exact field preservation, checksum and
duplicate rejection, and reviewed-byline/preprint/correction checks. The
relocation test exercises portable citation inputs, not an independently
installed or archived scientific pipeline.

## Original Export Limits

- The manually reviewed igraph citation is still in the service-reference
  supplement, not this CSL collection. This collection is the explicit
  selection, not an inferred full dependency inventory.
- The original OrthoHMM preprint Crossref record omits an author suffix;
  assembly preserves it unchanged. The previously documented source review
  and final journal rendering remain open.
- This artifact does not validate scientific claims, historical dataset
  releases, data redistribution rights or full transitive provenance.
- No native workflow, benchmark input, score, confidence interval or frozen
  scientific executor changed. No journal submission or external archive
  deposit was made.
