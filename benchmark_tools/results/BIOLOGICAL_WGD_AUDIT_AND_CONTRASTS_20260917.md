# WGD Application: Source Audit and Statistical Freeze

Supplement to `BIOLOGICAL_WGD_APPLICATION_PROTOCOL_20260917.md`, not a
replacement. Frozen after input mapping, before application inference or
inspection of any method's application outcomes. Input mapping is recorded in
`biological_wgd_inputs_20260917.json`; all240 source pairs remain retained.

## Chronology and Reference Scope

Scannell et al. (2011), DOI10.1534/g3.111.000273, explicitly place the WGD
in the ancestry of Saccharomyces sensu stricto (printed page19, Duplicate gene
losses). Their Figure1A marks its ancestral branch; Figure3 compares the focal
species. The abstract includes S.cerevisiae, S.mikatae, S.kudriavzevii and
S.bayanus var.uvarum; page12 explains their use of S.bayanus for the latter.
Thus the four-species common ancestor is post-WGD: retained WGD copies can be
ancestral paralogs at this root. This supports the chronological premise, not
every individual experimental pair's duplication assignment.

Primary source inspected: [author-hosted full article](https://dunham.gs.washington.edu/sensustricto.pdf).
The PMC viewer returned a browser check; the accessible author PDF was used.
No absolute divergence dates or strain equivalence are inferred here.

The local YGOBv7 README independently labels the four species and states that
pillar copy positions are arbitrary, without A/B-track or syntenic meaning.
Its SHA256 is162fb37cf0f7a50b44fae47f0d2c2f9b126b728660b65d9bee8d3175af8317f6.
Pillars support homology, not cross-species copy-specific orthology. Source
sequence/annotation errors and conflicting experimental assignments remain
possible. A missing/OFF protein is an input limitation, not a tool error.
The experimental classifications establish biological relevance, not orthology
truth. Development exposure and unresolved YGOB redistribution rights remain.

## Fixed Analysis Populations

Descriptive reporting retains all240pairs, including every exclusion. Among
239input-eligible pairs, report distinct groups, merged groups and incomplete
assignment. The231shared-unambiguous-pillar pairs define the fixed population
for all paired statistical contrasts below, independent of tool outcomes.
Also report each endpoint by High, Low and Sparse source class, descriptively
without additional intervals or hypothesis tests. Do not recode Sparse as Low.

For coverage only, exclude a pair if its reference pillar has no available
non-S.cerevisiae input member; determine and report this count before outcomes.
The remaining population is identical for all methods. Do not exclude pairs
because a method has missing predictions. A native run that fails admission
has unavailable endpoints, not fabricated zero scores. Preserve its failure
and do not replace it with a different configuration.

Input-only verification at freeze:231pairs occupy231distinct reference pillars;
none has a zero coverage denominator. Therefore all three contrast endpoints
use the same231pairs and231resampling units in this snapshot.

## Exact Paired Contrasts

Four oriented contrasts, each evaluated for the three endpoints below:

1. OrthoHMM satellite_v2 minus OrthoHMM high_sensitivity.
2. OrthoHMM satellite_v2 minus full OrthoFinder3.1.5.
3. OrthoHMM satellite_v2 minus SonicParanoid.
4. OrthoHMM high_sensitivity minus full OrthoFinder3.1.5.

The OrthoFinder MCL checkpoint is descriptive only. No comparisons are added
after results are inspected. No claim of a pure phylogeny effect follows from
contrast1: these frozen configurations differ in more than reconciliation.

Endpoint definitions for an admitted method:

- `separation_rate`: fraction of eligible pairs with both anchors explicitly
  assigned to distinct native groups. A missing anchor is zero, not a split.
- `supported_separation_rate`: fraction with distinct assigned anchor groups
  and at least one non-S.cerevisiae reference-pillar member in each group.
  Each group's support is checked independently. Missing anchors give zero.
- `mean_non_scer_coverage`: mean per-pair fraction of available non-S.cerevisiae
  reference-pillar genes in the union of the assigned anchor groups. A missing
  anchor contributes an empty group; both missing give zero. Denominators use
  reference/input membership, never a method's recovered set. Coverage can be
  high for a merged group and is not by itself a success criterion.

Groups mean root HOGs for the phylogenetic methods and native OGs for the
sequence methods, as in the original protocol. No synthetic singleton group
is created to turn an unassigned gene into a successful assignment. Duplicate
membership or unknown IDs fail output admission rather than selecting a group
arbitrarily. Explicit native singleton groups remain valid groups.

## Resampling and Multiplicity

Use the shared reference pillar as the resampling unit: all cohort pairs in
one pillar move together. Sort pillar IDs lexicographically and canonical ORF
pairs within each pillar. Draw the original number of pillars with replacement
20,000times using NumPy Generator(PCG64(20260920)). Share draw indices across
methods and endpoints. Recompute pair-weighted means from all pairs in each
sampled pillar, including multiplicities; do not average pillar means instead.
For coverage, remove only the outcome-independent zero-denominator pairs from
each replicate. Report any undefined replicates without replacing or redrawing
them; an endpoint with undefined replicates has no reported interval pending
an explicit statistical amendment, not a silently conditional interval.

Report point differences in percentage points. Use linear-interpolated
percentile intervals at0.025/0.975 and at0.05/(2*12),1-0.05/(2*12).
The12-comparison Bonferroni family remains12 even if a failed method makes
some comparisons unavailable. Adjusted percentile intervals are approximate
exploratory uncertainty, not guaranteed finite-sample coverage. Report the
number of pillars, pairs, ties and positive/negative per-pair differences.
Correlations beyond pillars and development exposure limit interpretation.
An interval spanning zero establishes neither equivalence nor superiority.

## Required Descriptive Guards

For every reference-eligible pair, retain assigned anchor-group sizes, native
group IDs, non-S.cerevisiae support per anchor, coverage numerator/denominator,
and the number of distinct native groups intersecting the reference pillar.
Report known foreign-pillar members and input members absent from the curated
reference separately within the union of assigned anchor groups, alongside
union size. Do not count unmapped members as verified contamination. Counts
are diagnostics, not orthology precision/recall, and get no extra inferential
tests. Include merged, unsupported, fragmented and incompletely assigned cases.

Keep all six previously hash-selected examples. YLR284C/YOR180C has conflicting
reference pillars: show that limitation rather than replacing the example.
The complete240-row supplement is mandatory. No composite score is optimized.

## Remaining Execution Gates

The chronology and exact contrast inventory are now specified. Before launch,
freeze commands, runtime identity and native-output admission for the four
methods using established publication settings on the original host. Recheck
the fixed coverage denominator population against the pinned input manifest. This
document does not itself authorize an unspecified native command. Do not use
the dedicated DGX while its scientific timing panel is active. Biological
outcomes, uncertainty, examples and figures remain outstanding.
