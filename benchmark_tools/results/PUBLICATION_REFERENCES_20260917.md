# Verified Publication References

Initial literature audit, checked 17 September 2026 local time. This is a
partial bibliography, not a completed citation or licensing audit. Method
papers document concepts; exact evaluated versions, commands, reference
snapshots and checksums are established by this project's run manifests.
No published benchmark numbers are substituted for the project's results.

## Benchmarks and Biological Evidence

1. Altenhoff AM et al. (2016). Standardized benchmarking in the quest for
   orthologs. Nature Methods 13:425-430.
   [DOI:10.1038/nmeth.3830](https://doi.org/10.1038/nmeth.3830).
   Publisher metadata and abstract checked. Supports the QfO benchmarking
   framework and its multiple tasks/trade-offs, not our selected six-metric
   mean or exact 2020 reference snapshot. See the local
   [scoring environment](qfo_assessment_environment_20260917.json).

2. Emms DM, Kelly S (2020). Benchmarking Orthogroup Inference Accuracy:
   Revisiting Orthobench. Genome Biology and Evolution 12:2258-2266.
   [DOI:10.1093/gbe/evaa211](https://doi.org/10.1093/gbe/evaa211).
   [Primary full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC7738749/)
   identifies the revised 70-RefOG benchmark and its supplied scoring suite.
   The [authors' repository](https://github.com/davidemms/Open_Orthobench)
   was also checked. Our retained exclusions and weighting still require
   the [local scoring protocol](ORTHOBENCH_UNCERTAINTY_PROTOCOL_20260916.md).

3. Byrne KP, Wolfe KH (2005). The Yeast Gene Order Browser: Combining
   curated homology and syntenic context reveals gene fate in polyploid
   species. Genome Research 15:1456-1461.
   [DOI:10.1101/gr.3672305](https://doi.org/10.1101/gr.3672305).
   Metadata verified on the [author's publication page](https://www.kevinbyrne.org/pubs/?field=gen);
   primary article text identifies curated sequence/synteny-based homology.
   This original paper describes seven species, not our later retained
   16-species subset. It does not establish family-disjoint validation or
   correct ancestral-copy labels for all current pillars. See the
   [frozen evaluation](YGOB_FROZEN_INTERPRETATION_20260916.md).

4. Manni M, Berkeley MR, Seppey M, Simao FA, Zdobnov EM (2021).
   BUSCO Update: Novel and Streamlined Workflows along with Broader and
   Deeper Phylogenetic Coverage for Scoring of Eukaryotic, Prokaryotic, and
   Viral Genomes. Molecular Biology and Evolution 38:4647-4654.
   [DOI:10.1093/molbev/msab199](https://doi.org/10.1093/molbev/msab199).
   [Author preprint](https://arxiv.org/abs/2106.11799) and
   [official documentation](https://busco.ezlab.org/busco_userguide) checked;
   publisher PDF metadata also identifies the journal citation. BUSCO
   addresses completeness using conserved ortholog expectations. Our
   Three Kingdoms pair score is a project evaluation over retained BUSCO
   labels, not a standard BUSCO completeness score or proteome-wide test.

5. Kuzmin E et al. (2020). Exploring whole-genome duplicate gene retention
   with complex genetic interaction analysis. Science 368:eaaz5667.
   [DOI:10.1126/science.aaz5667](https://doi.org/10.1126/science.aaz5667).
   [Primary full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC7539174/)
   and bibliographic record checked. Provides experimental evidence for
   240 duplicate pairs, not cross-species copy-specific orthology truth.
   The [local application protocol](BIOLOGICAL_WGD_APPLICATION_PROTOCOL_20260917.md)
   identifies the retained data deposit, table and prospective selection.

## OrthoFinder v3

6. Emms DM, Liu Y, Belcher L, Holmes J, Kelly S (2026). OrthoFinder:
   improved phylogenetic orthology inference with enhanced accuracy and
   scalability. Nature Methods 23:1327-1333.
   [DOI:10.1038/s41592-026-03126-6](https://doi.org/10.1038/s41592-026-03126-6).
   Publisher article, author list, publication date and change history
   checked. Relevant v3 method citation; it does not independently verify
   that our executable was version3.1.5, its child-tool paths, or our
   conversion semantics. Those remain run-provenance questions.

7. Emms DM, Liu Y, Belcher L, Holmes J, Kelly S (2026). Author Correction:
   OrthoFinder: improved phylogenetic orthology inference with enhanced
   accuracy and scalability. Nature Methods, published26August2026.
   [DOI:10.1038/s41592-026-03238-z](https://doi.org/10.1038/s41592-026-03238-z).
   Read the publisher correction: an incorrect Figure2 was replaced in
   the HTML/PDF article. No volume/page was supplied in the inspected
   correction citation, so none is invented here. Do not infer from this
   notice that our executed results changed or require replacement.

## Other Methods and Method Lineage

8. Li L, Stoeckert CJ Jr, Roos DS (2003). OrthoMCL: Identification of
   Ortholog Groups for Eukaryotic Genomes. Genome Research 13:2178-2189.
   [DOI:10.1101/gr.1224503](https://doi.org/10.1101/gr.1224503).
   [Primary full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC403725/)
   checked. Describes MCL-based grouping, not evidence that pre-MCL graph
   edges are final predictions. Our version1.4 runs and BLAST failure
   handling require the retained execution and conversion audits.

9. Cosentino S, Sriswasdi S, Iwasaki W (2024). SonicParanoid2: fast,
   accurate, and comprehensive orthology inference with machine learning
   and language models. Genome Biology 25:195.
   [DOI:10.1186/s13059-024-03298-4](https://doi.org/10.1186/s13059-024-03298-4).
   [Publisher article](https://link.springer.com/article/10.1186/s13059-024-03298-4)
   checked. The paper distinguishes species-pair relations and multispecies
   groups and describes multiple execution modes. It does not prove that
   all described components/settings were used in our2.0.9 runs. No paper
   speedup or accuracy estimate is transferred to our comparison.

10. Klemm P, Stadler PF, Lechner M (2023). Proteinortho6:
    pseudo-reciprocal best alignment heuristic for graph-based detection
    of (co-)orthologs. Frontiers in Bioinformatics 3:1322477.
    [DOI:10.3389/fbinf.2023.1322477](https://doi.org/10.3389/fbinf.2023.1322477).
    [Publisher article](https://www.frontiersin.org/journals/bioinformatics/articles/10.3389/fbinf.2023.1322477/full)
    checked. Appropriate major-version reference for retained6.3.6, with
    graph/search/clustering distinctions. Native post-clustering graph
    selection and actual command settings still require local provenance.

11. Majidian S, Nevers Y, Yazdizadeh Kharrazi A, Warwick Vesztrocy A,
    Pascarelli S, Moi D, Glover N, Altenhoff AM, Dessimoz C (2025).
    Orthology inference at scale with FastOMA. Nature Methods22:269-272.
    [DOI:10.1038/s41592-024-02552-8](https://doi.org/10.1038/s41592-024-02552-8).
    [Publisher PDF](https://www.nature.com/articles/s41592-024-02552-8.pdf)
    and [primary full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC11810774/)
    checked. Publication year is2025 despite2024 in the DOI. Our0.3.5
    results used a supplied OrthoFinder tree; this paper does not turn
    that run into independent tree inference or prove matching resources.

12. Emms DM, Kelly S (2015). OrthoFinder: solving fundamental biases in
    whole genome comparisons dramatically improves orthogroup inference
    accuracy. Genome Biology16:157.
    [DOI:10.1186/s13059-015-0721-2](https://doi.org/10.1186/s13059-015-0721-2).
    [Publisher article](https://link.springer.com/article/10.1186/s13059-015-0721-2)
    checked. Documents score normalization and orthogroup inference;
    use alongside the2019/2026 papers, not as a description of every
    operation or default in the current3.1.5 full pipeline.

13. Emms DM, Kelly S (2019). OrthoFinder: phylogenetic orthology
    inference for comparative genomics. Genome Biology20:238.
    [DOI:10.1186/s13059-019-1832-y](https://doi.org/10.1186/s13059-019-1832-y).
    [Publisher article](https://link.springer.com/article/10.1186/s13059-019-1832-y)
    checked. Documents the phylogenetic extension and distinct output
    levels. The sequence-only checkpoint and full phylogenetic outputs
    are not interchangeable merely because both yield gene groups.

14. Steenwyk JL, Buida TJ III, Rokas A, King N (2024). OrthoHMM:
    Improved Inference of Ortholog Groups using Hidden Markov Models.
    bioRxiv preprint.
    [DOI:10.1101/2024.12.07.627370](https://doi.org/10.1101/2024.12.07.627370).
    Metadata verified from the [author's publication list](https://jlsteenwyk.com/publications.html)
    and this repository's README/docs citation. Direct bioRxiv full-text
    retrieval failed during this check; no full-text verification is
    claimed. Preserve preprint status. This lineage citation does not
    document every built-in search or phylogenetic refinement change in
    the current frozen implementation; code and execution audits do.

## Outstanding Review

Comparator method metadata and the cited OrthoHMM preprint lineage are now
covered above. The [eight-dependency supplement](PUBLICATION_DEPENDENCY_REFERENCES_20260918.md)
adds simulator, sequence search, alignment, MCL and tree-method references.
The [CSL export](PUBLICATION_CITATION_EXPORT_20260918.md) now preserves complete
deposited author lists for these original 14 references. The
[numeric supplement](PUBLICATION_NUMERIC_REFERENCES_20260918.md) covers
Leiden/CPM and selected analysis dependencies; the
[service supplement](PUBLICATION_SERVICE_REFERENCES_20260918.md) covers
QfO2020/2022 and igraph, with explicit consortium-byline review flags.
HMM/profile implementation attribution, individual reference/functional
resources and journal-specific citation rendering remain incomplete.
Full-text review of the OrthoHMM preprint remains incomplete; its metadata
check is narrower than a review of its scientific claims or version history.

Website availability or an article's open-access label is not proof of
redistribution permission for every downloaded dataset, bundled executable,
annotation or derivative. Resource-level licensing, exact acquisition URLs,
retained-version identities and archival permissions remain separate work.
No new raw data were downloaded or redistributed during this citation check.
