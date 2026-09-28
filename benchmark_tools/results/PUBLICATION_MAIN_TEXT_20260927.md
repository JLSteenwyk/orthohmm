# OrthoHMM: HMM-Centered Group Inference With Phylogenetic Refinement

Condensed scientific draft, updated 28 September 2026. Not submission-ready. The
[extended manuscript](PUBLICATION_MANUSCRIPT_DRAFT_20260916.md) retains detailed
methods, historical analyses, citations and audit records. This main text
does not supersede frozen protocols or historical result manifests.

## Abstract

Orthology inference requires balancing homolog recovery against separation of
paralogs. We evaluated an HMM-centered pipeline with candidate-family expansion
and phylogenetic refinement against established orthology tools. QfO and
OrthoBench were the primary development-exposed benchmarks; a frozen YGOB
evaluation assessed transfer to additional taxa. Phylogenetic OrthoHMM improved
on its high-sensitivity configuration but did not consistently outperform full
OrthoFinder. On OrthoBench, their group-recovery F1 values were 74.11% and
72.74%, respectively. On corrected QfO, OrthoFinder had higher VGNC, SwissTrees
and TreeFam-A point estimates, whereas OrthoHMM had higher GO/EC similarity and
FAS. A matched-recall simulation control supported an initial HMM-search
contribution in a fixed graph procedure, not full-pipeline superiority.
Independent-family generalization, uncertainty for several QfO endpoints and
controlled comparative timing remain unresolved. The evidence supports a
bounded contribution and an explicit precision-recall trade-off, not universal
accuracy or efficiency claims.

## Methods

The retained configurations are high-sensitivity OrthoHMM and the satellite_v2
phylogenetic pipeline. The original OrthoHMM preprint [@orthohmm2024preprint]
describes its lineage, not all subsequent implementation changes. The prospective
method was frozen at `7f3a9e4`, with BLOSUM62, E-value threshold 1e-4, Leiden
CPM resolution 0.1 and seed 4 [@leiden2019; @cpm2011].
The phylogenetic configuration expands candidate families, infers gene and
species trees, and applies positive-paralogy pair inference. Historical runs
retain their actual revisions rather than inheriting this prospective pin.
The [method diagram](figures_publication_method_20260916/publication_method.pdf)
distinguishes initial search, profile refinement, candidate expansion and
phylogenetic inference.

Comparators were OrthoFinder 3.1.5 [@orthofinder2026; @orthofinder2026correction],
SonicParanoid 2.0.9 [@sonicparanoid2024], ProteinOrtho 6.3.6 [@proteinortho2023],
FastOMA 0.3.5 [@fastoma2025] and OrthoMCL 1.4 [@orthomcl2003]. These method
references do not replace run-specific executable provenance.
OrthoFinder's sequence-only MCL checkpoint was
a diagnostic output, not a separately finalized phylogenetic analysis.
FastOMA used a supplied OrthoFinder species tree. QfO inputs included native
ortholog pairs, native post-clustering relations or group-derived cross-species
pairs as appropriate; these are not interchangeable output semantics. The
[generated comparison](qfo_corrected_comparison_20260926_v7/scores.md)
reports each conversion and prediction count.

OrthoBench measures curated group recovery [@orthobench2020]. QfO
[@qfo2016; @qfo2020] reports GO and EC similarity,
VGNC, SwissTrees and TreeFam-A F1, and FAS separately. Their arithmetic mean
is a project-defined secondary summary, not an official QfO F1. Three Kingdoms
is supplementary BUSCO-reference recovery [@busco2021], not genome-wide orthology truth.
QfO and OrthoBench were repeatedly inspected during development. The
[YGOB protocol](YGOB_VALIDATION_PROTOCOL_20260916.md), using the curated
homology/synteny resource [@ygob2005], froze evaluation before
test-score inspection; overlap assessment limits its interpretation to
novel-taxon transfer rather than family-disjoint validation.

Uncertainty analyses resample declared reference units and recompute the
benchmark statistic within replicates. They retain prespecified multiplicity
adjustments, negative findings and unavailable contrasts. Curated families
are not automatically independent: shared evolutionary history and predictions
joining families can violate exchangeability. Dependent protein pairs are not
treated as independent observations to obtain narrower intervals.

## Results

### Accuracy Depends On The Endpoint

Phylogenetic OrthoHMM achieved 74.1061% OrthoBench F1 versus 72.7365% for full
OrthoFinder. A fresh installed full run reproduced all 59,770 OrthoHMM groups
exactly; this establishes reproducibility, not independent accuracy or a
controlled runtime comparison.
[Full-run verification](INTEGRATED_FULL_OB_RESULT_22337.md).

On corrected QfO, phylogenetic OrthoHMM versus full OrthoFinder scored 0.901690
versus 0.988546 on VGNC, 0.833513 versus 0.848413 on SwissTrees, and 0.614864
versus 0.791918 on TreeFam-A. OrthoHMM had higher GO similarity (0.490349
versus 0.469548), EC similarity (0.965650 versus 0.936130) and FAS (0.762993
versus 0.691422). These endpoints differ in reference scope and eligibility;
they do not form a single accuracy ranking.
[All eight methods](qfo_corrected_comparison_20260926_v7/scores.md).

For corrected SwissTrees, phylogenetic-minus-full-OrthoFinder F1 was -1.4900
percentage points, with a multiplicity-adjusted interval of [-8.7078, 7.2090].
Phylogenetic-minus-high-sensitivity OrthoHMM was +14.8015 points
[4.9791, 26.6015]. These conditional estimates use 18 development-exposed
families and a 24-endpoint adjustment; inclusion of zero is not equivalence.
Candidate expansion also changes between OrthoHMM configurations, so their
contrast does not isolate reconciliation.
[Protocol and contrast evidence](CORRECTED_SWISS_COMPARISON_RESULT_21987.md).

All 56 method-pair GO/EC comparisons had identical six-decimal scores on
shared scored pairs. Aggregate differences arose from different eligible
pair sets and their denominators. Restricting evaluation to intersections
would change the endpoint rather than resolve uncertainty. Historical raw
hashes were verified for all 16 input tables.
[Scored-pair audit](QFO_SCORED_PAIR_TRANSITIVE_BINDING_20260927.md).

### Component Evidence Is Bounded

The corrected QfO factorial supported candidate-expansion-by-reconciliation
interactions on SwissTrees F1, while all four profile-refinement F1 intervals
included zero. Reconciliation increased precision and reduced recall. These
results support conditional component effects, not proof that profile
expansion universally improves orthology inference.
[Complete factorial](QFO_CORRECTED_FACTORIAL_COMPLETE_20260919.md).

In 35 matched-recall simulation datasets spanning seven conditions and five
seeds, HMM-derived search evidence yielded mean downstream graph F1 of
83.6097%, versus 80.5180% for DIAMOND [@diamond2021]. The adjusted paired seed-block interval
for the +3.0917-point difference was [1.6782, 4.5313]. Profile expansion and
phylogeny were off. This is a development-exposed fixed-graph comparison,
not a comparison against OrthoFinder or a matched-effort result. Score
rankings and hit identities differ, so the design does not isolate a causal
mechanism. [Matched-recall control](MATCHED_GRAPH_RESULT_20260926.md).

### Transfer And Biological Recovery Reveal Trade-Offs

Frozen YGOB group-recovery F1 was 92.233654% for phylogenetic OrthoHMM and
92.318524% for full OrthoFinder. Their difference was -0.084870 points,
with adjusted interval [-0.622528, 0.445222]. OrthoHMM had higher precision
and lower recall. The evaluation projects predictions onto the reference
universe and assumes exchangeable pillars; it does not establish unrestricted
generalization or resolved pairwise orthology.
[YGOB results](YGOB_FROZEN_RESULTS_20260916.md).

A descriptive partition by the frozen overlap screen retained all original
false-positive allocations. In screen-negative pillars, phylogenetic OrthoHMM
versus full OrthoFinder had F1 of 57.93% versus 45.02%, precision of 55.92%
versus 30.74%, and recall of 60.10% versus 84.03%. In screen-positive pillars,
their F1 values were 93.37% versus 95.01%. The negative stratum contained
2,893 singleton pillars out of 3,298, compared with 2,017 out of 6,952 in the
positive stratum. This composition difference and the screen's inability to
exclude remote homology prevent independent-family or causal interpretations.
The [all-method stratum figure](figures_ygob_overlap_20260928/ygob_overlap_strata.pdf)
is descriptive, with no new uncertainty estimate or method tuning.

The prespecified whole-genome-duplicate application used experimental evidence
from Kuzmin and colleagues [@kuzmin2020]. All 240 experimental pairs were
retained, with 239 input-eligible and 231 shared-reference-pillar pairs.
On this development-exposed application, phylogenetic OrthoHMM separated
238 of 239 input-eligible pairs, versus 58 for high sensitivity and 236 for
full OrthoFinder. Among the 231 reference-eligible pairs, however, only 193
OrthoHMM separations retained at least one non-S. cerevisiae reference homolog
with each anchor, versus 56 for high sensitivity and 227 for full OrthoFinder.
Mean per-pair homolog coverage was 82.338%, 99.149% and 98.413%, respectively.
Coverage counts reference homologs in the union of the anchor groups and can
be high even when the anchors are merged; it is not orthology recall.
The [complete five-method comparison](BIOLOGICAL_WGD_RESULTS_20260917.md)
also retains SonicParanoid and the diagnostic OrthoFinder MCL checkpoint.

Phylogenetic OrthoHMM minus full OrthoFinder had a supported-separation
difference of -14.719 percentage points, with adjusted interval
[-21.645, -8.225], and a coverage difference of -16.075 points
[-19.755, -12.496]. These exploratory percentile intervals use 20,000 paired
pillar resamples and a 12-endpoint adjustment. Each eligible pair occupies
a distinct pillar; these conditional estimates do not establish independent
generalization. Supported separation does not establish
cross-species copy-specific orthology. In the
[six prespecified case traces](BIOLOGICAL_WGD_CASE_TRACE_20260917.md), five
focal homologs left their anchor groups during root-lineage reconstruction,
before satellite constraints. They remained in other output groups; no example
was replaced to improve the result.

A subsequent [fixed-tree diagnostic](WGD_FIXED_TREE_RULE_RESULTS_20260928.md)
evaluated all four existing root-duplication rules on the same seven candidate
families and six examples, with its protocol frozen before alternative outcomes
were examined. The supported-children and confidence rules reproduced the
species-overlap baseline partitions exactly. The mapped-event rule reduced
homolog coverage in four of five reference-eligible examples, without improving
supported separation. None recovered the five focal homologs into either anchor
group. Native pair predictions and confidence annotations remained unchanged.
This post hoc intervention rejects these three alternatives as a repair on the
fixed inputs, not as methods in general. It neither establishes topology error
nor identifies ancestral-copy truth, and no defaults were changed.

## Discussion And Limitations

The supported contribution is an HMM-centered alternative with measurable
component effects and explicit accuracy trade-offs. Initial HMM search has
bounded support in matched-recall simulations; additional profile refinement
has not demonstrated a general benefit. Neither simulation evidence nor an
OrthoBench point advantage establishes superiority over full OrthoFinder.

Original TreeFam-A family mappings and complete source trees remain unavailable.
Public archive recovery yielded historical Selectome subtrees, but all are
restricted to Euteleostomi. A [taxonomic coverage audit](TREEFAM_RECOVERED_SCOPE_20260928.md)
found that 55,933 of 79,320 retained reference relations (70.52%) have at least
one endpoint outside that clade. Those relations cannot be reconstructed from
these subtrees under a species-consistent mapping; the other 23,387 relations
are only potentially in scope, not demonstrated reconstructions. This is a
coverage exclusion, not an accuracy effect or a basis for family-level intervals.
VGNC dependence and
rare-error diagnostics do not justify the candidate confidence-interval
procedure. Valid paired uncertainty for GO/EC, FAS and the secondary mean is
also unfinished. These gaps must remain visible alongside point estimates.
The [corrected FAS sample audit](QFO_CORRECTED_FAS_SAMPLE_AUDIT_20260928.md)
reproduces all eight means, but scored fractions range from 0.0067% to 58.46%
of reported eligible pairs and every sample reuses proteins across pairs.
Native pair-level standard errors do not resolve this comparison uncertainty.
The completed simulations and tree perturbations do not cover arbitrary
evolutionary conditions, and novel-taxon YGOB testing retains homolog-family
overlap with development data.

Recovered OrthoMCL results retain sequence-specific BLAST failures. Failed
high-CPM experiments remain unavailable, not zero-scoring or successful runs.
Shared-host runs and the historical DGX panel are descriptive resource records;
the replacement controlled Threadripper timing panel has not run. No controlled speedup is
claimed. Versioned archival release, complete dependency/data redistribution
review, journal formatting and final visual review remain open.

## Reproducibility And Availability

The extended manuscript links frozen protocols, source revisions, seeds,
input/output manifests, conversion audits, executable workflows and generated
figures. Installed inference and independent readback have been validated on
full OrthoBench, but that is not proof of all-method cross-host portability.
Historical locks are retained as provenance; patched replacement environments
must be distinguished from the binaries used for original scores. The
[progress ledger](PUBLICATION_PROGRESS.md) records completed work and unmet
requirements. No submission-ready release or archival DOI is claimed.

## References

The bibliography attributes the methods and resources cited above. Full
execution evidence and additional dependency references remain in the extended
manuscript; this selected list is not complete software attribution.
