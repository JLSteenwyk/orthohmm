# Reference-Resource Citation Supplement

Checked 18 September 2026. Resource descriptions do not authenticate the
exact reference snapshot, scoring implementation, or redistribution rights.
This supplements, rather than replaces, retained input and scoring manifests.

## TreeFam

Li H et al. (2006). TreeFam: a curated database of phylogenetic trees of
animal gene families. Nucleic Acids Research 34:D572-D580.
[DOI:10.1093/nar/gkj118](https://doi.org/10.1093/nar/gkj118).
The [publisher article](https://academic.oup.com/nar/article/34/suppl_1/D572/1133584)
was inspected for bibliographic metadata and its distinction between
automatically generated TreeFam-B and manually curated TreeFam-A trees.
This paper describes release1.1, not the release7 inputs identified by the
2016 QfO paper. It supplies resource attribution, not provenance for our
2020 pooled reference or proof that family identities can be reconstructed.
The [source-retrieval investigation](TREEFAM_SOURCE_RETRIEVAL_20260918.md)
still records the missing original trees and mapping.

## SwissTree

SIB Swiss Institute of Bioinformatics. SwissTree. Resource entry,
[ExPASy](https://www.expasy.org/resources/swisstree), accessed18September2026.
The official description identifies a curated reference-gene-phylogeny
resource intended for assessing phylogenomic databases. Its stated aim of
100 gold-standard phylogenies is not the size of our retained benchmark.
Our SwissTrees analysis uses18 families; the local count audits and scoring
reference establish that subset. This is a website citation, not an invented
journal citation or evidence that all reference families are independent.

## Feature Architecture Similarity

Dosch J, Bergmann H, Tran V, Ebersberger I (2023). FAS: assessing the
similarity between proteins using multi-layered feature architectures.
Bioinformatics 39(5):btad226.
[DOI:10.1093/bioinformatics/btad226](https://doi.org/10.1093/bioinformatics/btad226).
Metadata checked against the [bibliographic record](https://pubmed.ncbi.nlm.nih.gov/37084276/);
the [publisher's indexed article](https://academic.oup.com/bioinformatics/article/39/5/btad226/7135831)
describes architecture comparison and resolution of overlapping annotations.
Direct full-text retrieval failed during this check; a complete full-text
review is not claimed. The [authors' software documentation](https://bionf.github.io/FAS/)
independently describes the feature-architecture score and overlap handling.
FAS similarity is not pairwise orthology F1. This2023 citation does not prove
that the QfO2020 feature archive used that implementation or default settings.
Our [retained-sample audit](QFO_FAS_SAMPLE_AUDIT_20260917.md) remains the
evidence for reproduced aggregate arithmetic, not independent feature
annotation validation or dependency-aware uncertainty.

## Remaining Work

The [selected Crossref export](publication_resource_references_20260918.csl.json)
now contains TreeFam and FAS with
[source provenance](publication_resource_citation_provenance_20260918.json).
Checksum-verified cached replay reproduces its bytes exactly. Crossref's
TreeFam record contains onlyH.Li; this omission is not treated as the true
article byline. The separate
[reviewed export](publication_resource_bylines_20260918.csl.json) supplies
all15authors from DOI-checked article XML, matching the publisher byline.
[Correction provenance](publication_resource_byline_provenance_20260918.json)
records the raw/export/XML/parser hashes. Only that author field changed;
the FAS record is unchanged.26export/parser tests pass.

Reproduce the byline transformation from the retained inputs with the existing
`benchmark_tools.correct_qfo_service_bylines.byline(xml_bytes, doi)` function:
read the raw CSL, replace its `Li2006TreeFam` author field with the returned
list from `PMC1347480.xml`, and serialize using
`json.dumps(records, indent=2, ensure_ascii=True) + "\n"`. Verify the recorded
input hashes first; the expected output hash is in the correction provenance.
This article XML is not the missing TreeFam-A tree/mapping archive.

VGNC, GO/EC and other annotation-resource attribution, complete machine-readable
export of the SwissTree website citation, journal-specific formatting, and exact-version
and rights review remain separate requirements. No benchmark inputs, scores,
confidence intervals or method settings changed during this citation review.
