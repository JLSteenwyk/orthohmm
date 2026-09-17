# Prospective WGD Paralog Application

Freeze this protocol before mapping the experimental cohort to any tool outputs.
The source study, table schema and source class counts have been inspected;
family-level OrthoHMM/comparator outcomes for this application have not. This
is an application of the frozen method, not another independent generalization
test. Saccharomyces sequences/families overlap development resources and were
excluded from the earlier novel-taxon YGOB experiment for that reason. Do not
retune the method or claim family independence from this application.

## Independent Cohort

Use all240 distinct paralog pairs in TableS7 of Kuzmin et al., Science2020,
DOI10.1126/science.aaz5667, PMID32586993. The source experimental classification
has47high-fraction,114low-fraction and79sparse-interaction pairs. Preserve all
three strata; do not select only experimentally well-connected or famous genes.
The experimental measurements motivate biological relevance and define strata,
not orthology truth. High/low labels are source measurements, not universal
functional equivalence/divergence labels.

Sources: [primary study](https://pubmed.ncbi.nlm.nih.gov/32586993/),
[deposited data](https://doi.org/10.5061/dryad.g79cnp5m9),
[Zenodo mirror](https://zenodo.org/records/3975054).
Use the August2020 deposit represented by Zenodo3975054. TableS7.xlsx:
39,665bytes, depositedMD5cefe40a24fd97cf135d8f3827dcd0d65, SHA256
bf502168d9c87956edf6c764650e45a8c60a8605e08ea93fd7ef57a49313f7fe.
Retained API metadata SHA256
a113657e4a5b914fb2c06ddaa3fc35c941b995c7e3fc1c3db00d7d6a4c3222e8.
The deposit declares CC0; the article's license is separate. Do not redistribute
YGOB sequences until its independent redistribution terms are resolved.

## Input and Reference Gates

Use complete available YGOBv7 proteomes for S.cerevisiae, S.mikatae,
S.kudriavzevii and S.bayanus var.uvarum. Their common ancestor postdates the WGD,
so retained WGD copies are ancestral paralogs for this input scope. Verify that
chronology against primary phylogenetic evidence before inference/admission.
No cohort-only sequence subset is supplied to inference. Snapshot complete
input IDs, byte hashes, species assignments, and duplicates/missing sequences.

Map experimental systematic ORF IDs exactly, without sequence-based rescue or
manual outcome-driven remapping. Preserve all240cohort rows with explicit
mapping/exclusion reasons. For family-coverage diagnostics require both anchors
to map unambiguously to the same YGOB homology pillar. Conflicts, absent anchors,
ambiguous IDs or disputed WGD assignments remain separate unresolved categories.

YGOB's README explicitly says positions1/2 are arbitrary and are not ancestral
A/B tracks. Never infer cross-species copy-specific orthology from column order.
The homology pillars are curated using sequence and genomic context, not wholly
sequence-independent truth. Any later copy-specific positive benchmark requires
separately curated syntenic/phylogenetic evidence and a new prospective protocol.

## Frozen Comparisons

Run frozen OrthoHMM high_sensitivity and satellite_v2, full OrthoFinder3.1.5 and
SonicParanoid on the same four complete proteomes. Retain the OrthoFinder MCL
checkpoint as a diagnostic, not a separately run sequence-only3.1.5 pipeline.
Use established publication configurations and admit source/runtime/commands
and complete outputs. Record failures and coverage. Run on the original host
without competing on the dedicated DGX timing machine; do not claim controlled
runtime comparisons for this application.

For phylogenetic methods evaluate root-level HOG membership at this four-species
scope; for sequence methods use native orthogroups. Explicitly record the output
semantics and extraction, especially OrthoFinder's root N0 HOG table. Native
cross-species pair tables alone cannot answer a same-species separation test.

## Endpoints and Safeguards

Primary descriptive endpoint: among mapped experimental pairs, report both
anchors assigned to distinct groups, both merged into one group, and incomplete
assignment separately. A missing/unassigned anchor is not a successful split.
Always also report the full240-pair denominator and all input/reference exclusions.

Singleton splitting is insufficient biological evidence. Report a stricter
supported-separation endpoint: anchors in distinct groups and each group also
contains at least one non-S.cerevisiae member of the anchors' curated homology
pillar. This is homolog support, not proof of correct copy-specific orthology.
Alongside it report non-S.cerevisiae pillar coverage across the union of both
anchor groups, foreign-pillar contamination and unmapped predicted members.
Keep these separate endpoints; do not optimize a composite score or call them
orthology precision/recall. Report both clean splits and overfragmentation,
merges and unsupported recovery. Retained root-HOG groups need not reproduce
the older pre-WGD homology pillar as a single group.

Use all source strata, with sparse measurements retained as unknown functional
class rather than zero redundancy. If intervals are reported, pair methods on
the same eligible units, recompute each statistic, and group experimental pairs
sharing a reference pillar into one resampling unit. Use20,000PCG64 draws,
seed20260920, nominal95% intervals and a declared correction over all selected
method/endpoint contrasts; freeze the exact contrast inventory before outcomes.
Intervals are exploratory on development-exposed data, not an independent
superiority or causal-function test.

## Examples and Completion

Before outcomes, select two examples per source stratum by ascending SHA256 of
the canonical ORF pair joined with a tab (six total). Do not replace examples
because a tool fails or the plot is uninteresting. Show excluded examples with
their reason. Report the complete240-row supplement so failures are visible.
Literature-supported cases may be discussed additionally, clearly labeled as
post hoc and never substituted for the prospective examples.

Completion requires a frozen cohort, input/reference audit, verified native
outputs, all endpoints and exclusions, uncertainty/contrast inventory where
used, a complete case table and figures. This protocol alone does not establish
biological usefulness. Do not claim correct cross-species copy assignment or
general orthology improvement from a same-species split plus homolog coverage.
