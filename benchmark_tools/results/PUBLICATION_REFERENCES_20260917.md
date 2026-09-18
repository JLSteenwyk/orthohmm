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

## Outstanding Review

Still verify the OrthoHMM publication lineage; OrthoFinder's preceding
methods papers; OrthoMCL, SonicParanoid2, Proteinortho6 and FastOMA; simulator,
HMM/search, clustering and phylogeny dependencies; QfO service updates and
individual reference/functional resources. Complete author lists and a
journal-formatted citation export remain to be produced.

Website availability or an article's open-access label is not proof of
redistribution permission for every downloaded dataset, bundled executable,
annotation or derivative. Resource-level licensing, exact acquisition URLs,
retained-version identities and archival permissions remain separate work.
No new raw data were downloaded or redistributed during this citation check.
