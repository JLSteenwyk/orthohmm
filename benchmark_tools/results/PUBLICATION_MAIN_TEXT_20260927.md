# OrthoHMM: HMM-Centered Group Inference With Phylogenetic Refinement

Condensed scientific draft, 27 September 2026. Not submission-ready. The
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
phylogenetic pipeline. The prospective method was frozen at `7f3a9e4`, with
BLOSUM62, E-value threshold 1e-4, Leiden CPM resolution 0.1 and seed 4.
The phylogenetic configuration expands candidate families, infers gene and
species trees, and applies positive-paralogy pair inference. Historical runs
retain their actual revisions rather than inheriting this prospective pin.
The [method diagram](figures_publication_method_20260916/publication_method.pdf)
distinguishes initial search, profile refinement, candidate expansion and
phylogenetic inference.

Comparators were OrthoFinder 3.1.5, SonicParanoid 2.0.9, ProteinOrtho 6.3.6,
FastOMA 0.3.5 and OrthoMCL 1.4. OrthoFinder's sequence-only MCL checkpoint was
a diagnostic output, not a separately finalized phylogenetic analysis.
FastOMA used a supplied OrthoFinder species tree. QfO inputs included native
ortholog pairs, native post-clustering relations or group-derived cross-species
pairs as appropriate; these are not interchangeable output semantics. The
[generated comparison](qfo_corrected_comparison_20260926_v7/scores.md)
reports each conversion and prediction count.

OrthoBench measures curated group recovery. QfO reports GO and EC similarity,
VGNC, SwissTrees and TreeFam-A F1, and FAS separately. Their arithmetic mean
is a project-defined secondary summary, not an official QfO F1. Three Kingdoms
is supplementary BUSCO-reference recovery, not genome-wide orthology truth.
QfO and OrthoBench were repeatedly inspected during development. The
[YGOB protocol](YGOB_VALIDATION_PROTOCOL_20260916.md) froze evaluation before
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
83.6097%, versus 80.5180% for DIAMOND. The adjusted paired seed-block interval
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

The prespecified whole-genome-duplicate application assessed paralog separation
together with retention of reference homologs. Phylogenetic OrthoHMM's greater
separation than high sensitivity did not remove its homolog-recovery losses.
Full OrthoFinder remained an important negative comparison. Raw separation
alone is insufficient because splitting off anchors can increase separation
without useful ortholog recovery. The complete application, controls and
subsequent traces remain in the extended manuscript rather than selecting
only favorable examples.

## Discussion And Limitations

The supported contribution is an HMM-centered alternative with measurable
component effects and explicit accuracy trade-offs. Initial HMM search has
bounded support in matched-recall simulations; additional profile refinement
has not demonstrated a general benefit. Neither simulation evidence nor an
OrthoBench point advantage establishes superiority over full OrthoFinder.

Original TreeFam-A family mappings remain unavailable. VGNC dependence and
rare-error diagnostics do not justify the candidate confidence-interval
procedure. Valid paired uncertainty for GO/EC, FAS and the secondary mean is
also unfinished. These gaps must remain visible alongside point estimates.
The completed simulations and tree perturbations do not cover arbitrary
evolutionary conditions, and novel-taxon YGOB testing retains homolog-family
overlap with development data.

Recovered OrthoMCL results retain sequence-specific BLAST failures. Failed
high-CPM experiments remain unavailable, not zero-scoring or successful runs.
Shared-host runs and the historical DGX panel are descriptive resource records;
the replacement dedicated timing panel has not run. No controlled speedup is
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
