# Clustering And Analysis References

This supplement adds five explicitly selected references, distinct from
the existing comparator/benchmark and simulator/external-tool bibliographies.
DOI, title, author and publication-year fields were retrieved through the
Crossref CSL metadata service. Publisher/project sources confirm their roles:

| Reference | Role in this work | Source |
| --- | --- | --- |
| Traag, Waltman and van Eck (2019), From Louvain to Leiden: guaranteeing well-connected communities | Leiden optimization used by frozen clustering | [Paper](https://doi.org/10.1038/s41598-019-41695-z) |
| Traag, Van Dooren and Nesterov (2011), Narrow scope for resolution-limit-free community detection | CPM objective used by `CPMVertexPartition`; not a separate inferred-orthology guarantee | [Paper](https://doi.org/10.1103/PhysRevE.84.016114), [package documentation](https://leidenalg.readthedocs.io/en/stable/advanced.html) |
| Harris et al. (2020), Array programming with NumPy | Numerical arrays and the paired-bootstrap implementation | [Paper](https://doi.org/10.1038/s41586-020-2649-2) |
| Cock et al. (2009), Biopython: freely available Python tools for computational molecular biology and bioinformatics | FASTA and phylogenetic-tree parsing in validation tools | [Project citation listing](https://biopython.org/wiki/Publications), [paper](https://doi.org/10.1093/bioinformatics/btp163) |
| Hunter (2007), Matplotlib: A 2D Graphics Environment | Publication plots | [Project citation guidance](https://matplotlib.org/stable/project/citing.html), [paper](https://doi.org/10.1109/MCSE.2007.55) |

The frozen `7f3a9e4` implementation's `orthohmm/externals.py` imports igraph
and leidenalg and invokes `find_partition` with `CPMVertexPartition`.
`bootstrap_qfo_factorial.py` uses NumPy; `validate_factorial_native.py` uses
Biopython `Phylo` and `SeqIO`; `plot_qfo_factorial.py` uses Matplotlib.
These source observations establish why the references are included, not
which package version was executed in every historical run. Exact versions
remain bound by each run's runtime and environment manifests.

## Artifacts And Checks

- [Selection](publication_numeric_citation_selection_20260918.json)
- [CSL-JSON](publication_numeric_references_20260918.csl.json), SHA-256
  `d9d21dc677fba634463bb4f8b52f6c0323e046a38514a18abf6c2b8caf8c0ee2`
- [Retrieval provenance](publication_numeric_citation_provenance_20260918.json), SHA-256
  `e63fc1137f0027efb4fee76b73c551b8cdcd13f91f76711944c98212b8a2d509`

Raw metadata responses are retained locally in
`benchmarks/work/publication_numeric_crossref_20260918/`. A network-free
re-export using the saved provenance reproduces the committed CSL bytes
exactly. Article identifiers5233 and016114 remain CSL `number` fields,
not invented page ranges. Biopython's metadata date is its online publication
date, not a claim that its print issue appeared in March.

This is a bibliographic supplement, not full-text review, version provenance,
license clearance or proof of numerical correctness. No HMMER citation is
used here to imply that the built-in search kernel executes HMMER. igraph is
now covered by the [service/library supplement](PUBLICATION_SERVICE_REFERENCES_20260918.md).
Additional reference resources, remaining dependencies and journal-specific
rendering still need their own coverage in the final bibliography.
