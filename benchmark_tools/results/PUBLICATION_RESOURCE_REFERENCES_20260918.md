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

## VGNC

Jones TEM, Yates B, Braschi B, Gray K, Tweedie S, Seal RL, Bruford EA
(2023). The VGNC: expanding standardized vertebrate gene nomenclature.
Genome Biology24:115.
[DOI:10.1186/s13059-023-02957-2](https://doi.org/10.1186/s13059-023-02957-2).
The [primary article](https://link.springer.com/article/10.1186/s13059-023-02957-2)
describes consensus orthology predictions from Ensembl Compara, NCBI Gene,
PANTHER and OMA informing nomenclature assignment, with manual review in
specified cases. This is relevant reference-design context, not proof that
every pair in the retained QfO snapshot used that exact workflow or an
estimate of comparator-specific circularity. It does not establish independent
experimental orthology truth. The2023resource paper does not authenticate
the2020benchmark snapshot.

## Gene Ontology

Ashburner M et al. (2000). Gene Ontology: tool for the unification of biology.
Nature Genetics25:25-29.
[DOI:10.1038/75556](https://doi.org/10.1038/75556).

The Gene Ontology Consortium (2026 issue; online18December2025). The Gene
Ontology knowledgebase in2026. Nucleic Acids Research54(D1):D1779-D1792.
[DOI:10.1093/nar/gkaf1292](https://doi.org/10.1093/nar/gkaf1292).
The [current GO citation policy](https://geneontology.org/docs/go-citation-policy/)
requests the original and current resource papers, together with the actual
data release/version. These citations do not imply that our QfO2020 scorer
used current annotations, and are not citations for its similarity formula.
Exact historical annotation and ontology identities remain a separate
provenance requirement.

The [publisher's correction](https://academic.oup.com/nar/article/54/12/gkag678/8721936)
addresses consortium author Daiqing Chen's name. The retrieved Crossref
record already contains that spelling. It includes a consortium literal
plus135individual authors; a full author-by-author review is not claimed.
Its `issued` field is18December2025, whereas the volume issue date is
6January2026. The raw CSL preserves the online date rather than silently
changing it to match the title year. Journal-specific rendering remains
separate.

## Annotation Citation Export

The [three-record annotation export](publication_annotation_references_20260918.csl.json)
and [provenance](publication_annotation_citation_provenance_20260918.json)
retain selected DOI metadata, source-response hashes and retrieval times.
Cached replay reproduces the CSL bytes exactly;16exporter tests pass. No
benchmark data, scores or method configuration changed.

```sh
python benchmark_tools/export_publication_citations.py \
  --manifest benchmark_tools/results/publication_annotation_citation_selection_20260918.json \
  --raw-directory benchmarks/work/publication_annotation_crossref_20260918 \
  --output NEW_CSL.json --provenance NEW_PROVENANCE.json \
  --cached-provenance benchmark_tools/results/publication_annotation_citation_provenance_20260918.json
```

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

EC and other annotation-resource attribution, complete machine-readable
export of the SwissTree website citation, journal-specific formatting, and exact-version
and rights review remain separate requirements. No benchmark inputs, scores,
confidence intervals or method settings changed during this citation review.
