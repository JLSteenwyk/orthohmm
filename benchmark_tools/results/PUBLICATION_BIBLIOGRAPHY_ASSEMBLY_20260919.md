# Consolidated Selected Bibliography

## Metadata Review Rendering

The [rendered bibliography](publication_bibliography_review_20260919/bibliography.html)
contains all 37 current v3 records. Pandoc 3.1.3 citeproc rendered the original
project `publication-review.csl` style, not an assumed target journal style.
The retained Pandoc AST and manifest record the commands, tool/source/input
hashes, output hashes, empty stderr and per-entry field-presence checks.

Reproduce from the repository root with a new output directory:

```sh
python benchmark_tools/render_publication_bibliography.py \
  --bibliography benchmark_tools/results/publication_bibliography_20260919_v3.csl.json \
  --style benchmark_tools/publication-review.csl --output NEW_REVIEW_DIRECTORY
python -m pytest -q tests/unit/test_render_publication_bibliography.py
```

All 19 renderer tests pass, including actual 37-record rendering after input
relocation, unchanged CSL input bytes, output checksums, inventory mismatch
rejection, hidden link-target exclusion, and detection of missing titles,
identifiers, consortium names, family names, suffixes, preprint labels and
dates. The development render initially demoted family-name particles;
the style now explicitly retains them adjacent to the family name. The
failed development output remains outside the publication results in
`benchmarks/work/publication_bibliography_review_initial_20260919/`.

Style SHA-256:
`8fcc6447a351b2b9b3612af121f9c259eb3bf32a608d0b9548d7d5580e4ef4a8`.
Rendered HTML SHA-256:
`8e58768b467b15c01ebce304fe5a041d52c2eb30bbb9bbdb8d1f629eb19cad8b`.

The audit checks inventory and substring visibility, not full citation
semantics or visual layout. Publisher, volume, issue, pagination and given-name
typography are rendered but not independently validated by these checks.
This is not journal formatting, complete manuscript citation coverage, rights
clearance, or validation of scientific results. No benchmark input or method
changed. The earlier export-history sections below retain their original
then-outstanding limitations.

## Current Author-Reviewed Export

`publication_bibliography_20260919_v3.csl.json` retains all 37 records and
changes exactly one field from v2: the second OrthoHMM author's CSL suffix
is `III`. The [author's publication list](https://jlsteenwyk.com/publications.html)
explicitly includes this suffix for DOI 10.1101/2024.12.07.627370. The
[bioRxiv API](https://api.biorxiv.org/details/biorxiv/10.1101/2024.12.07.627370)
omits it in both deposited versions, as does the retained Crossref record.
This is an author-source-supported rendering correction, not a claim that
the publisher supplied the suffix. The article page returned HTTP 403;
no access restriction was bypassed and no full-text byline verification is
claimed. Article date, author order, title and preprint status are unchanged.

The correction is explicit in
`publication_orthohmm_suffix_provenance_20260919.json`, including before/after
author fields and SHA-256-bound source snapshots. Raw Crossref and all
earlier exports remain unchanged. The locally retained author page and API
response have SHA-256 values respectively
`d713be57a6457ba281f9a275f290df253ab095692de61c9aba4c58431b1d5a25`
and `fbbf4f22c8b2f98bada38e48ffdd7e1eb1efb5864edb6c03ffad4737bcf0d933`.
These source webpages/API responses are not committed.

```sh
python benchmark_tools/correct_orthohmm_citation_suffix.py \
  --raw benchmark_tools/results/publication_references_20260918.csl.json \
  --author-page benchmarks/work/steenwyk_publications_20260919.html \
  --api benchmarks/work/orthohmm_biorxiv_metadata_20260919.json \
  --output NEW_REVIEWED_REFERENCES.csl.json --provenance NEW_SUFFIX_PROVENANCE.json
python benchmark_tools/assemble_publication_citations.py \
  --manifest benchmark_tools/results/publication_bibliography_selection_20260919_v3.json \
  --output NEW_BIBLIOGRAPHY.csl.json --provenance NEW_PROVENANCE.json
```

The first command reproduces the reviewed 14-record source; the second uses
the retained, checksum-bound copy selected in v3. Final bibliography SHA-256:
`353ff5508b647d8b49bdd5ca64b49f5ec870a4000f2dd0de73a942e234f68c4a`.
Forty-one focused tests pass, including source-change rejection, output
overwrite refusal, exact one-field correction and relocated v3 assembly.
The suffix review is now resolved at the stated author-source evidence level;
journal rendering, full bibliography coverage and rights review remain open.

## igraph Export (v2)

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
