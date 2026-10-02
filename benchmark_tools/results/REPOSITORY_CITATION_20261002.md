# Repository Citation Metadata

Add root [CITATION.cff](../../CITATION.cff) and a README link. This prepares
machine-readable citation metadata, not a release or archival deposition.
The software author follows setup.py. The preferred citation follows the
[reviewed reference](publication_references_reviewed_20260919.csl.json) and
[author-suffix provenance](publication_orthohmm_suffix_provenance_20260919.json):
the original four-author 2024 bioRxiv preprint, not the current benchmarking
manuscript. CFF `unpublished` is used for its preprint status.

The [author publication list](https://jlsteenwyk.com/publications.html)
supports the Buida III suffix omitted by the deposited Crossref metadata.
The [Crossref DOI record](https://api.crossref.org/works/10.1101/2024.12.07.627370)
was read successfully through curl; its author order, title, DOI and year
agree with the retained review. Web-tool access to that API failed; do not
claim the web-tool fetch succeeded. No full-text or new scientific claims
were validated. GitHub documents root CFF files and preferred citations in
its [citation-file guide](https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/about-citation-files).

## Executed Validation

Use a new isolated environment, not a frozen scientific runtime:

```bash
python -m venv citation-validation
citation-validation/bin/python -m pip install cffconvert==2.0.0
citation-validation/bin/cffconvert --validate --infile CITATION.cff
```

Actual private environment:
`benchmarks/work/publication_citation_20261002/venv`, Python 3.12.3.
Validation exits zero against the bundled CFF 1.2.0 schema. Structured YAML
parsing and that schema's Draft7 validator also pass eight checks: schema,
title, ordered names/suffix, DOI/URL, year, explicit preprint/software-DOI
separation, absent invented release fields, and setup.py software author.
Two negative schema cases reject missing title and unsupported `preprint`
type. These are focused metadata checks, not a scientific test-suite run.

Local BibTeX and APA-like conversion commands exit zero, but cffconvert 2.0.0
renders only the top-level software entry and ignores preferred-citation.
Inspect those outputs as software rendering only. The preferred fields were
checked through parsed metadata, not GitHub's rendered citation UI.

The [validation receipt](repository_citation_validation_20261002.json)
records package versions, schema/input/output identities and all checks.
Receipt: 3,148 bytes; SHA256
`b9be395fe382c94a22a28761f2db1a511508f69c96fd9fd98a2b2de3c027634d`.
An exploratory API inspection used a nonexistent Citation.as_bibtex method
and raised AttributeError; it did not run validation or alter metadata.
A patch newline typo was corrected before validation.

## Release Boundary

Do not invent a software DOI, released version/date, or benchmarking-paper
author list. Readers must report exact revision and inference settings.
The preprint DOI belongs only to preferred-citation, not software metadata.
The existing annotated v0.5.0 tag resolves to historical commit
`1ea3d2fae0d47edc94dfee65153905a373393adb`, not frozen scientific revision
`7f3a9e40dd7e79f842cc2c11fb8b548f9a802806`; neither tag nor revision changes.
This does not assert that the historical tag is the publication release.

No inference, timing observation, archive rebuild, manuscript rerender,
DGX access, service change, shared-environment upgrade, release upload or
deposition executes. Controlled resource evidence remains outstanding;
publication readiness is not established.
