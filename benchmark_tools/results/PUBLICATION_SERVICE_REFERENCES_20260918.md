# QfO Service And igraph Citation Coverage

## Selected References

- Altenhoff AM et al. (2020). The Quest for Orthologs benchmark service and
  consensus calls in 2020. Nucleic Acids Research48(W1):W538-W545.
  [DOI](https://doi.org/10.1093/nar/gkaa308),
  [primary article](https://pmc.ncbi.nlm.nih.gov/articles/PMC7319555/).
  Supports service/reference-proteome history and multi-test evaluation,
  not our exact corrected input archive or project-defined six-metric mean.
- Nevers Y et al. (2022). The Quest for Orthologs orthology benchmark service
  in 2022. Nucleic Acids Research50(W1):W623-W632.
  [DOI](https://doi.org/10.1093/nar/gkac330),
  [primary article](https://pmc.ncbi.nlm.nih.gov/articles/PMC9252809/).
  Describes the 2020 reference-proteome update and VGNC benchmark, as well as
  endpoint-specific sensitivity/specificity trade-offs. It does not establish
  our results or a single universally best method. The service currently
  [recommends citing this update](https://orthology.benchmarkservice.org/proxy/).
- Csardi G, Nepusz T (2006). The igraph software package for complex network
  research. InterJournal, Complex Systems,1695.
  [Official Python-interface citation guidance](https://python.igraph.org/en/0.11.6/).
  This is the graph-library citation, distinct from the Leiden optimizer and
  CPM objective. No DOI for the 2006 article is asserted here; the C-library
  software archive DOI is not silently substituted for an article DOI.
  The documentation version used as citation evidence does not establish the
  igraph version executed by any retained run.

## Metadata And Reproducibility

The two QfO records were exported with the existing Crossref exporter:
`publication_service_citation_selection_20260918.json`,
`publication_service_references_20260918.csl.json`, and
`publication_service_citation_provenance_20260918.json`.
The CSL SHA-256 is
`8c844fbef9a2694da284a2e98e22baf78a5b865465bcabd7768b321e877fc53a`.
Raw responses remain in `benchmarks/work/publication_service_crossref_20260918`.
A cached, network-free re-export reproduces the CSL byte-for-byte;16 exporter
tests pass. Deposited online dates are distinct from print-issue dates.

The igraph guidance HTML is retained locally at
`benchmarks/work/igraph_citation_guidance_20260918.html` (12143bytes), SHA-256
`51683b3ae3f119aea2e2a12fb8200fc87679e5cbc6b076035ad34d8b04f39c7a`.
Its recommended bibliographic fields were checked directly. The webpage is
not redistributed in the source repository. The citation above is manually
transcribed guidance, not a Crossref record or full-text review.

## Reviewed Byline Export

The source discrepancy is now resolved in a separate CSL export,
`publication_service_bylines_20260918.csl.json`, with a complete before/after
author record in `publication_service_byline_provenance_20260918.json`.
The raw Crossref export above remains unchanged. All non-author CSL fields
are identical. Use the corrected export for these two publication citations.

The article XML snapshots were downloaded from Europe PMC's fullTextXML
endpoint for [PMC7319555](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7319555/fullTextXML)
and [PMC9252809](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9252809/fullTextXML)
into `benchmarks/work/publication_service_crossref_20260918/`.
SHA-256 values are respectively
`50162a498e935d807943160e706d3ed3f6f6f5769313277e9de4eea0fe4074ec`
and `d6b9ae310a6fd7bf617d6011e4649313dbe1992b4df5dd76af1ff57231f50255`.
They are locally retained evidence, not redistributed article text.

`correct_qfo_service_bylines.py` parses only the top-level article-meta
contributor group, verifies DOI and snapshot checksums, preserves contributor
order and personal-name spelling, and does not deduplicate people by name.
The 2020 list remains23 entries, with the consortium's given/family encoding
converted to a CSL literal. The 2022 list changes from69 expanded Crossref
entries to31 top-level byline entries. Its XML combines OpenEBench and QfO
in one collective label; that label is retained verbatim apart from whitespace,
not arbitrarily split into two invented contributor records.

26 focused tests pass (10 byline-parser and16 existing citation-export tests).
Fresh-path offline replay reproduces the corrected CSL byte-for-byte; a direct
comparison confirms every non-author field remains unchanged. This resolves
the author-list rendering issue, not full journal-style bibliography formatting.

### Original Discrepancy

Preserve the Crossref export as deposited evidence, not a final journal
bibliography: the2020 consortium name is split into given/family fields;
the2022 record combines collective names and expands consortium members,
including repeated personal names and spelling variants. Publisher byline
and consortium membership must be distinguished in a separately documented
rendering correction; do not silently deduplicate the raw author list.
The checked2022 article byline identifies OpenEBench and QfO as collective
authors, with consortium membership listed separately in the article.

Full journal-style rendering, remaining functional/reference resources and
third-party rights remain separate requirements. These references do not
attribute the built-in HMM kernel to HMMER or certify executed versions.
