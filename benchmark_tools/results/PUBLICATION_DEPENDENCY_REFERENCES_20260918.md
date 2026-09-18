# Dependency Literature Audit

This supplements `PUBLICATION_REFERENCES_20260917.md`. Eight selected
DOI records were exported with the tested citation exporter. Exact local
versions, executable resolution and actual invocation remain execution
provenance questions; the literature is not evidence that a particular
binary or optional module ran. No published speedup is imported into our
runtime comparison.

| Dependency | Primary literature | Scope in this project |
| --- | --- | --- |
| Zombi | [Davin et al., Bioinformatics 36:1286-1288](https://academic.oup.com/bioinformatics/article/36/4/1286/5578480), DOI 10.1093/bioinformatics/btz710 | Species/genome histories and gene trees; pinned revision and seed adapter are specified in the simulation protocol. |
| Pyvolve | [Spielman and Wilke (2015)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0139047), DOI 10.1371/journal.pone.0139047 | Sequence simulation along gene trees; retained version 1.1.0 and per-family seeds are local execution choices. |
| MAFFT | [Katoh and Standley (2013)](https://doi.org/10.1093/molbev/mst010) | Multiple alignment method reference, not proof of a particular strategy or thread count. |
| FastTree | [Price, Dehal and Arkin (2010)](https://pmc.ncbi.nlm.nih.gov/articles/PMC2835736/), DOI 10.1371/journal.pone.0009490 | Approximate maximum-likelihood tree inference; distinguish retained 2.1.11 and 2.2.0 executables. |
| DIAMOND | [Buchfink, Reuter and Drost (2021)](https://doi.org/10.1038/s41592-021-01101-x) | Protein search and sensitivity modes; literature does not establish matched sensitivity with OrthoHMM. |
| FAMSA | [Deorowicz, Debudaj-Grabysz and Gudys (2016)](https://www.nature.com/articles/srep33964) | Alignment method lineage; this paper alone does not document all changes in bundled 2.2.3. |
| FastME | [Lefort, Desper and Gascuel (2015), author reference page](https://www.atgc-montpellier.fr/fastme-papers-contact-sup-mat/), DOI 10.1093/molbev/msv150 | Distance-based phylogeny method; installation/resolution does not establish invocation by every comparator run. |
| MCL | [Enright, Van Dongen and Ouzounis (2002)](https://academic.oup.com/nar/article/30/7/1575/2376029), DOI 10.1093/nar/30.7.1575 | Protein-family clustering method; not equivalent to resolved ortholog-pair inference. |

Publisher or author primary pages were inspected for all eight entries.
The Zombi publisher page gives online publication 30 September 2019 and
issue date February 2020. Crossref's `issued` is 2019, which the unmodified
CSL export preserves; the final journal style must consistently handle the
2020 issue citation. The deposited title says "dead linages"; do not silently
alter an article title to match the differently spelled repository heading.
Full deposited author names, including accents, are retained in CSL-JSON.

The Zombi paper describes a multilevel simulator and validation of event
distributions. Our panel does not exercise all its features: transfer and
species extinction were disabled; the variable-length panel varies length
between families, not through within-family indels. Simulator citation does
not replace the local event-truth conversion tests or validate a new model.
See `PUBLICATION_SIMULATION_PROTOCOL_20260916.md`,
`PUBLICATION_VARIABLE_LENGTH_PROTOCOL_20260916.md`, and
`run_zombi_seeded.py`. Native ARM companion-tool resolution and limitations
are recorded in `DGX_SCALING_MIGRATION_20260917.md`.

## Export And Remaining Work

Selection: `publication_dependency_citation_selection_20260918.json`.
Output: `publication_dependency_references_20260918.csl.json`.
Provenance: `publication_dependency_citation_provenance_20260918.json`.
Raw responses: `benchmarks/work/publication_dependency_citations_20260918/`.
Use `export_publication_citations.py --cached-provenance` with these paths
and fresh output destinations for checksum-verified offline replay.

The earlier 14-reference export remains unchanged. This is an additional
eight-record set, not a replacement or duplicate acquisition of that set.
HMM implementation/profile libraries, Leiden, statistical/plotting
dependencies, functional annotation resources, and resource-specific reuse
licenses still need coverage. FAMSA version-specific follow-up and final
journal citation rendering remain open. Do not call this the complete
publication bibliography or a complete software license inventory.
