# OrthoHMM: HMM-Centered Group Inference With Phylogenetic Refinement

Condensed scientific draft, native-evidence revision of 6 October 2026. Not submission-ready.
This revision preserves the
[third 4 October text](PUBLICATION_MAIN_TEXT_20261004_v3.md), its 35-page assembled
review and the rc4 archive. Fresh native ablations, the two-cell QfO figure
and functional-pair composition below are later additions, not retroactive
members of those artifacts. Native inference and remaining publication work
are incomplete; the progress ledger and actual scheduler state govern execution.
The prior revision's chronology remains recorded below.

This evidence-integrated working revision preserves the
[second 4 October text](PUBLICATION_MAIN_TEXT_20261004_v2.md), its 34-page review
and the rc3 archive. The newly integrated inventory, candidate-cap mechanism
and retained stage-cost addenda are not retroactively part of those artifacts.
It also preserves the
[first 4 October text](PUBLICATION_MAIN_TEXT_20261004.md) and its rendered reviews;
it does not relabel their historical resource statements or change the frozen
method, benchmark scores or endpoints. The
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
Family-disjoint generalization and uncertainty for several QfO endpoints
remain unresolved. The completed shared-host matched-resource panel retains
all 27 attempts, with 24 eligible measurements; it does not establish isolated
comparative efficiency. The evidence supports a
bounded contribution and an explicit precision-recall trade-off, not universal
accuracy or efficiency claims. Fresh native ablations retain neutral
downstream-profile effects and a SwissTrees precision-recall trade-off whose
adjusted F1 interval includes zero. Higher native GO/EC means reflect scored-pair
membership differences, not changed serialized similarities on shared pairs.
Synthetic null scoring revealed
composition-dependent behavior of the approximate significance filter.
Heterogeneous-length simulations favored full OrthoFinder; supplied generating
trees did not establish an accuracy gain, and stronger topology perturbations
reduced F1 and recall in several conditions. A fixed-candidate generating
gene-tree diagnostic localized most divergent-condition recall loss upstream
of reconciliation; residual event-history traces identified satellite-constraint
losses and duplication evidence absent from retained candidates.

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

Initial search uses banded local maximum-path scoring with match, insert and
delete states, integer BLOSUM62 emissions and fixed transition costs.
Raw scores receive an approximate significance filter before geometric-mean
length normalization. MSA-derived profiles learn match emissions, while
retaining uniform transition costs; strict expansion requires both a member-score
threshold and an initial sequence-supported anchor. This is not demonstrated
HMMER/phmmer equivalence or calibrated assignment confidence. The
[numerical specification](FROZEN_HMM_SCORING_SPECIFICATION_20260930.md)
records the recurrence, constants, profile construction and acceptance rules
at the full frozen source revision. No scoring/default change accompanies it.

A [prespecified null-score protocol](FROZEN_NULL_SCORE_PROTOCOL_20261002.md)
examined that significance filter without fitting parameters. Three residue
compositions and lengths 50/150/400 each used ten seeds and 1,000 independently
drawn query-target pairs per seed. Each pair had a one-target database and was
scored with full and width-64 bands: 90,000 independent pairs, not 180,000
independent observations. Five fixed cutoffs yielded 90 tail endpoints with
Bonferroni-adjusted exact binomial intervals. These intervals describe the
synthetic generator, not reference-family or orthology uncertainty.

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
QfO and OrthoBench were repeatedly inspected during development. A retained
family-evidence inventory scanned 1,785 immutable JSON result blobs at Git
revision `f145fff1` and 114 separately hash-frozen local-only reports, including
historical parameter sweeps. The exact 70 OrthoBench/18 retained SwissTrees
identities were canonicalized; score blocks, named/list mentions and unresolved
ordinal/synthetic labels remained separate. Every block retained its file,
JSON pointer, count-vector digest and recorded partition. Reports were not
counted as independent native experiments or evidence of causal tuning.
[Inventory protocol](DEVELOPMENT_FAMILY_INVENTORY_PROTOCOL_20261004.md).
The [YGOB protocol](YGOB_VALIDATION_PROTOCOL_20260916.md), using the curated
homology/synteny resource [@ygob2005], froze evaluation before
test-score inspection; overlap assessment limits its interpretation to
novel-taxon transfer rather than family-disjoint validation.

Uncertainty analyses resample declared reference units and recompute the
benchmark statistic within replicates. They retain prespecified multiplicity
adjustments, negative findings and unavailable contrasts. Curated families
are not automatically independent: shared evolutionary history and predictions
joining families can violate exchangeability. Dependent protein pairs are not
treated as independent observations to obtain narrower intervals.

The prespecified QfO parameter neighborhood retained seven arms: the frozen
control and six one-at-a-time changes. CPM resolutions were 0.08/0.12 versus
0.1, candidate minimum normalization 0.024/0.036 versus 0.03, and minimum
margin 1.2/1.8 versus 1.5. SwissTrees F1 was recomputed from macro precision
and recall in each of 100,000 shared family-bootstrap draws, not averaged
across family F1 values. The seed was 20260925; Bonferroni adjustment retained
all 18 planned endpoints across six contrasts and three metrics. This is
development-exposed local sensitivity, not independent method selection.
[Frozen parameter protocol](QFO_PARAMETER_NEIGHBORHOOD_PROTOCOL_20260919.md).

Two evolutionary simulation panels used Zombi and Pyvolve
[@zombi2019online; @pyvolve2015], ten seeds and seven conditions each: baseline,
divergence, duplication/loss turnover, their combination, missing proteins,
uneven sampling and a taxon-count control. Event-derived cross-species truth
was not replaced by ancestral-family membership. Fixed 300-residue sequences
and heterogeneous ancestral-family lengths of 100-500 residues were analyzed
separately; descendants retain family length, without within-family indels.
Native admission, including finite comparator graph weights, preceded scoring.
Within each panel, 20,000 paired successful-seed resamples recomputed mean
seed-level metrics, with a 14-contrast adjustment for primary F1. Failures were
not imputed; available-case means from different seed sets are not paired effects.
[Simulation protocols and results](SIMULATION_VARIABLE_NATIVE_INTERPRETATION_20260916.md).

On the 70 heterogeneous-length datasets, phylogenetic OrthoHMM and full
OrthoFinder also used supplied generating trees and deterministic rooted
nearest-neighbor interchanges at clade distances two and four. Generating
trees are oracle diagnostics, not achievable end-to-end inference or guaranteed
accuracy upper bounds. All 126 F1/precision/recall endpoints were exploratory,
with 20,000 paired-seed resamples, seed 20260918 and Bonferroni adjustment.
[Tree-control protocol](SIMULATION_TREE_CONTROL_PROTOCOL_20260917.md).

The subsequent fixed-candidate gene-tree diagnostic retained all 70 OrthoHMM
generating-species-tree cells. Only inferred, single-ancestral-family candidates
received generating gene trees, either at their generating root or rerooted by
the frozen minimum-duplication/loss procedure. Candidates, species trees,
reconciliation and satellite-constraint policy were fixed; native pairs had to
reproduce exactly before intervention. Bypasses and mixed-ancestry candidates
retained native predictions. Descriptive ten-seed means and counts were reported
without population intervals or independent-confirmation claims. The oracle
changes branch lengths/support as well as topology/root, and is not an accuracy
upper bound or a topology-only causal control.
[Gene-tree diagnostic protocol](SIMULATION_GENE_TREE_ORACLE_PROTOCOL_20261004.md).

A retained upstream trace followed all 163,527 true pairs through significant
initial HMM hits, direct final graph edges, undirected connectivity, final
candidates and native pairs. The saved graph precedes refinement/expansion;
it is neither initial search evidence nor ortholog output. Missing significant
hits do not distinguish prefilter rejection, score filtering, caps or ranking.
Gene-to-seed membership inside merged candidates was not reconstructed from
incomplete sidecars. Independent igraph/SciPy component and native-row checks
covered all 70 cells.
[Upstream trace protocol](SIMULATION_UPSTREAM_TRACE_PROTOCOL_20261004.md).

Finally, a complete post hoc cohort selected every generating-root candidate
with within-candidate FP/FN after screening all 10,125 candidates. All
cross-species pairs in those candidates, including correctly classified pairs,
were traced against original simulator event tables and reconciled XML.
Induced ancestor clades, species overlap, raw pair calls, root groups and native
satellite-constraint filtering had to reproduce original truth and oracle counts.
An independent XML/Biopython reader recalculated fixed-rule and constraint
decisions without importing the trace worker. No defaults or endpoints changed.
[Residual trace protocol](SIMULATION_ORACLE_RESIDUAL_PROTOCOL_20261004.md).

A post hoc reproducibility diagnostic compared accepted candidate-merge traces
from the original OrthoBench satellite arm and three already completed native
repeats by round and named source/target memberships. All common cluster
indices and forward/reverse hit counts were checked. A controlled 19-gene
fixture used the unchanged frozen candidate engine, four attachments per
anchor per round and two rounds, varying only initial cluster order or one
satellite's scores by one floating-point unit. Rejected candidates and the
origin of historical numerical differences were not reconstructed.
[Candidate-cap diagnostic](CANDIDATE_TRACE_VARIATION_RESULT_20261004.md).

### Fresh Native Ablation Protocol

A separate frozen execution panel comprises six fresh OrthoBench identities
and seven fresh corrected-QfO identities. P denotes downstream profile expansion,
C candidate-family expansion and R reconciliation; initial sensitive HMM search
is on in every cell. The retained P1/C0/R0 configuration is a separately reused
historical baseline, not an additional fresh run. Native inference uses the
frozen method with fresh output directories and does not reuse phylogenetic
checkpoints. These identities are distinct from the completed 27-attempt scaling
panel. Terminal scheduler, scientific-output and resource-accounting reviews
remain separate gates; recovery of valid scientific outputs after a wrapper or
measurement failure does not repair the original timing result.
[Frozen native plan](native_factorial_receipt_amendment_20261004/plan.json).

Fresh RefOG or SwissTrees count vectors and metadata must match retained
family records exactly before their original paired intervals can be reused.
OrthoBench retains all 12 planned simple effects and the 36-endpoint adjustment,
with 20,000 original draws (seed 20260918). QfO retains all 14 planned contrasts
and the 42-endpoint adjustment, with 100,000 original draws (seed 20260922).
No new bootstrap draws, independent confirmation, family-disjoint validation
or total-HMM-versus-non-HMM comparison follows from this binding. Reconciliation
also changes the QfO output from group-clique relations to inferred ortholog pairs.

A read-only functional diagnostic binds six native raw GO/EC/FAS tables through
their original scientific admissions. It preserves original means and
denominators while comparing shared and exclusive scored pairs. Independent
CSV/Decimal parsing and SQLite joins check the actual tables without importing
the primary comparison parsers. These are composition and dependence diagnostics,
not confidence intervals, corrected endpoints or another annotation-scoring run.
[Diagnostic protocol](NATIVE_QFO_FUNCTIONAL_PAIR_PROTOCOL_20261006.md).

### Shared-Host Resource Measurement

The replacement resource panel uses the local x86 Threadripper host, bizon,
not the historical DGX/ARM environment. The frozen panel comprises
high-sensitivity OrthoHMM, satellite_v2 OrthoHMM with inferred phylogeny and
full OrthoFinder 3.1.5 on nested sets of 4, 8 and 12 proteomes, with three
rotated repeats per method and size (27 sequential attempts). The input counts
are 73,266, 165,168 and 251,378 proteins. Scientific settings, input bytes and
order remain frozen. Native tasks use the same 32-CPU affinity and 128-GiB
ceiling; each non-exclusive Slurm allocation reserves 64 slots to accommodate
measurement and reporting work. The selected affinity uses distinct physical
cores, but allocation does not reserve the host against unrelated processes.

Under the user-authorized 3 October amendment, background analyses do not
invalidate a run merely by competing for resources. Launch capacity, source
and input identity, process attribution, measurement continuity and failure
handling still apply. Process snapshots and pressure observations bracket the
native interval and retain background CPU demand, memory pressure and I/O
diagnostics. Snapshot coverage is not continuous isolation certification, and
CPU pressure can exceed nominal bounds because of non-atomic accounting;
pressure magnitudes are diagnostic rather than timing-eligibility thresholds.

The primary endpoints are native-command monotonic elapsed time, native
task-subtree CPU-stat change including wrapper work, and native-step lifetime
memory peak including the launcher. They exclude separate preparation,
conversion and scoring; the peak is not pure algorithm RSS. Terminal scheduler,
runtime, environment, resource and native-output evidence is reviewed before
the next identity. All attempts and failures are retained without fastest-run
selection, automatic retries or estimated overhead subtraction. Median/range
summaries require all three eligible repeats; observed ranges are not confidence
intervals. Shared-host distortion is unknown and potentially method dependent.
Historical DGX and earlier shared-host timings are not pooled with this panel.

Timing measurements were collected on a shared Threadripper while other analyses
were running. Competition for CPU, memory bandwidth and I/O may have affected
elapsed times, with an unknown and potentially tool-dependent impact. These
are observed shared-host timings, not estimates of isolated performance.
The 6 October authorization does not require the DGX or a quiet window;
safe capacity and valid accounting remain launch prerequisites.
[Execution amendment](PUBLICATION_SHARED_HOST_AMENDMENT_20261003.md),
[endpoint and continuation contract](THREADRIPPER_SHARED_CONTINUATION_20261003.md).

## Results

### Accuracy Depends On The Endpoint

Phylogenetic OrthoHMM achieved 74.1061% OrthoBench F1 versus 72.7365% for full
OrthoFinder. In the exploratory eight-method comparison, the F1 difference
was +1.370 percentage points (21-endpoint adjusted interval [-7.281, 12.166]).
The precision difference was +15.705 [1.175, 30.634] and recall difference
-13.151 [-26.947, -0.609]. This is a precision-recall tradeoff, not established
F1 superiority. The complete panel uses 100,000 paired RefOG draws and
remains conditional on family exchangeability and development exposure.
[All 21 contrasts](OB_COMPLETE_UNCERTAINTY_RESULT_20260928.md).
The [complete interval figure](figures_ob_complete_uncertainty_20260928/ob_complete_uncertainty.pdf)
shows F1, precision and recall on the same difference scale without ranking
methods by their observed effects.
A fresh installed full run reproduced all 59,770 OrthoHMM groups
exactly; this establishes reproducibility, not independent accuracy or a
controlled runtime comparison.
[Full-run verification](INTEGRATED_FULL_OB_RESULT_22337.md).
A second full run using a separately reconstructed base Python and fresh
inference/reader environments also reproduced all groups and all 70 family
score records. Native pair, confidence, reconciliation-event and hierarchy
files were byte-identical to the earlier run. This extends same-host
environment reproduction, not biological validation or complete archive
restoration. [Reconstructed-base verification](RECONSTRUCTED_FULL_OB_RESULT_22376.md).
A third full run used restored local execution assets and separately acquired
upstream inputs. Its independent scientific admission reproduced all 59,770
groups, all 70 family score records and four native TSVs exactly, without
reusing phylogenetic checkpoints. This validates local archive-to-results
reproduction, not cross-host/OS restoration, independent accuracy or controlled
timing. [Restored-archive verification](RESTORED_ARCHIVE_FULL_OB_RESULT_22377.md).

On corrected QfO, phylogenetic OrthoHMM versus full OrthoFinder scored 0.901690
versus 0.988546 on VGNC, 0.833513 versus 0.848413 on SwissTrees, and 0.614864
versus 0.791918 on TreeFam-A. OrthoHMM had higher GO similarity (0.490349
versus 0.469548), EC similarity (0.965650 versus 0.936130) and FAS (0.762993
versus 0.691422). These endpoints differ in reference scope and eligibility;
they do not form a single accuracy ranking.
[All eight methods](qfo_corrected_comparison_20260926_v7/scores.md).

An exploratory [VGNC deletion diagnostic](CORRECTED_VGNC_INFLUENCE_RESULT_20260928.md)
removed each of 16,844 reference blocks and its incident scored pairs in turn,
across all eight methods. Every comparator retained a negative F1 difference
against full OrthoFinder for every deletion. For phylogenetic OrthoHMM the
range was [-8.7353, -8.6729] percentage points. These are sensitivity ranges,
not confidence intervals: dependent deletions change the fixed scored table,
and single-block sign stability does not establish population or joint-deletion
robustness. Native inference and eligibility were not rerun.

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
The [complete pair-composition figure](figures_qfo_scored_pair_decomposition_20260930/qfo_pair_decomposition.pdf)
displays the shared-denominator and exclusive-pair terms for all seven
comparators against full OrthoFinder. Phylogenetic OrthoHMM's shared pairs
account for 90.12% of its GO scored set but 46.40% of OrthoFinder's, and
86.59% versus 58.00% for EC. The original rounded-mean differences, +2.080130
and +2.952029 score points, are the sums of much larger opposing terms.
These fractions describe eligible scored pairs, not proteome-wide coverage;
the decomposition is arithmetic, not causal attribution or paired uncertainty.

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

### Fresh Native OrthoBench Effects

All six fresh OrthoBench cells have admitted scientific scores over 70 RefOGs;
420 family records match the retained factorial exactly. Partition agreement
is not asserted for every non-reference gene. P0/C0/R0 retains failed wrapper
22427 with scientific recovery, not successful timing. The two omitted native
cells and their contrasts remain unavailable rather than imputed.

| Native Cell | Group-Recovery F1 (%) | Precision (%) | Recall (%) |
| --- | ---: | ---: | ---: |
| P0/C0/R0 | 69.7634 | 78.8686 | 62.5429 |
| P0/C0/R1 | 72.7050 | 87.4948 | 62.1922 |
| P0/C1/R0 | 66.6395 | 64.4883 | 68.9391 |
| P0/C1/R1 | 73.4023 | 81.4862 | 66.7776 |
| P1/C0/R1 | 73.3114 | 87.4586 | 63.1038 |
| P1/C1/R0 | 67.2074 | 64.7001 | 69.9169 |

[Native point estimates](native_factorial_progress_20261005_v6/report.json),
[conditional interval binding](native_factorial_uncertainty_binding_20261005_v4.json),
[independent family-record readback](native_factorial_uncertainty_readback_20261005_v2.json).

Six of 12 planned simple effects have matching native records. With C0/R1 fixed,
profile expansion changes F1 by +0.6064 percentage points (adjusted interval
[-1.2114, 4.5210]; family wins/ties/losses 4/61/5). With C1/R0 fixed, its effect
is +0.5680 [-1.1132, 4.2906] (6/59/5). Both intervals include zero. These
are downstream-profile effects with initial HMM search retained, not the total
HMM contribution or proof of general benefit.

With P0 fixed, reconciliation changes F1 by +2.9416 [0.2542, 8.1443] at C0
and +6.7628 [1.3762, 14.6116] at C1. Candidate expansion changes F1 by
-3.1239 [-10.6588, 3.2751] at R0 and +0.6973 [-2.4590, 4.8209] at R1,
while increasing recall and reducing precision in both cases. These conditional
development-exposed group-recovery contrasts retain family-exchangeability and
percentile-coverage limitations. They neither complete a fresh eight-cell
factorial nor establish superiority over full OrthoFinder. Group co-membership
is not resolved pairwise-ortholog truth.
[Interpretation and exact endpoints](NATIVE_FACTORIAL_UNCERTAINTY_BINDING_20261005_V4.md).

### Fresh Native QfO Trade-Off And Pair Composition

Two of seven fresh corrected-QfO cells have admitted accuracy; five remain
unavailable in this dated snapshot. Both shown cells have P0/C0, with initial
HMM search on. They are not the selected-default all-tool comparison above.
The recovered R1 scientific result retains original failed 22437 timing,
null resources and timing ineligibility. Successful conversion, assessment and
scientific admission do not turn that failed timing into a valid measurement.

| Endpoint | Statistic | P0/C0/R0 | P0/C0/R1 |
| --- | --- | ---: | ---: |
| VGNC | F1 | 0.66683353 | 0.89818458 |
| SwissTrees | F1 | 0.68918396 | 0.78957434 |
| TreeFam-A | F1 | 0.60540356 | 0.60250783 |
| GO | Similarity | 0.47211905 | 0.49025964 |
| EC | Similarity | 0.93211406 | 0.96770222 |
| FAS | Sample mean | 0.77494473 | 0.78500682 |

[Admitted native endpoint snapshot](native_qfo_scientific_scores_20261006_v1/report.json).
The corresponding project-defined secondary means are 0.69009982 and
0.75553924, not official QfO F1 or a single accuracy ranking. Submitted relations
decrease from 9,009,082 group-derived clique pairs to 5,113,820 inferred pairs;
the fractions of 984,137 input accessions appearing in predictions are 55.5258%
and 55.1107%. Mapping loss is zero in both conversions. These are prediction
coverage and output-semantics observations, not accuracy statistics.

Only the R effect at P0/C0 has exact matched SwissTrees family records in both
fresh native cells; the other 13 planned contrasts remain unavailable. Reusing
the original draws and all 42 adjusted endpoints gives F1 +10.0390 percentage
points with adjusted interval [-4.6418, 24.4199], precision +30.5212
[14.3252, 46.5473] and recall -6.5334 [-18.0413, -0.0436].
F1 wins/ties/losses are 14/1/3. The F1 interval includes zero; precision rises
while recall falls. Original count-based arithmetic differs slightly from
native serialized endpoints and remains separately reported, not substituted.
The 18 families are development-exposed, and exchangeability and percentile
coverage are conditional assumptions.
[Guarded interval binding](native_qfo_swiss_uncertainty_binding_22449_20261006.json),
[independent arithmetic readback](recovered_native_qfo_swiss_readback_22449_20261006.json).

**Native QfO Figure.** The [reviewed four-panel figure](native_qfo_p0c0_figure_20261006_v1/native_qfo_p0c0.pdf)
separates the three F1 endpoints from GO/EC/FAS, shows SwissTrees precision and
recall, and retains the zero-crossing adjusted F1 effect. The small TreeFam-A
marker overlap does not conceal its decrease in the table above. Values come
directly from the admitted snapshot and guarded binding; exact original-Python
3.10 replay is separate from Python 3.12 rendering. This is a two-cell ablation,
not an initial-HMM-off control, a complete fresh factorial or a timing figure.
[Figure source/readback/visual review](NATIVE_QFO_FIGURE_RESULT_20261006.md).

Every R-on GO scored pair (78,607) occurs in the R-off set (145,142), and every
R-on EC scored pair (116,929) occurs in R-off (186,098). All common serialized
six-decimal scores agree exactly. The 66,535 GO and 69,169 EC R-off-only pairs
have lower means, 0.4506870677 and 0.8719528955. Thus, at retained precision,
the higher native functional means reflect scored-pair membership and
denominator composition, not changed similarities on common pairs or proof of
correct ortholog selection. Intersection-only means would change the endpoint.

FAS realizes 38,205 R-off and 252,451 R-on sampled pairs, with only 1,007 shared
pairs: 2.6358% of R-off and 0.3989% of R-on. All shared serialized scores agree;
the original sample-mean difference remains +0.0100620873, not replaced by a
shared-only difference of zero. Unseeded method-specific sample mixtures,
missing scores and repeated proteins prevent this overlap from supplying an
admitted paired confidence interval. Independent SQLite readback checks all
817,432 rows in the six bound tables. Neither membership arithmetic nor native
pair-IID standard errors resolve other-endpoint comparison uncertainty.
[Functional composition result](NATIVE_QFO_FUNCTIONAL_PAIR_RESULT_20261006.md),
[independent SQL readback](native_qfo_functional_pair_sql_readback_20261006.json).
No new score, default, independent-validation claim or timing repair accompanies
this diagnostic. GO/EC/FAS, VGNC, TreeFam-A and secondary-mean paired uncertainty
remains unfinished.

### Synthetic Null Tails Depend On Composition

At the frozen `E < 1e-4` cutoff, width-64 scoring passed 0/1/0 of 10,000
ordinary-background pairs at lengths 50/150/400. These sparse counts do not
establish rare-tail calibration. Half-background-plus-half-glutamine sampling
passed 89.83%, 100% and 100%; its adjusted intervals were [88.7484%, 90.8453%]
and [99.9181%, 100%] for each longer length. All
[90 endpoints and their figure](figures_frozen_null_scores_20261002_v2/frozen_null_scores.pdf)
are retained. A hit here is only a forced pair passing the approximate
significance gate, not a predicted ortholog or an observed pipeline false
positive. No prefilter or biological inference ran. The strong synthetic
composition is not a measured real-proteome prevalence. The
[audit](FROZEN_NULL_SCORE_RESULT_20261002.md) reproduces all endpoint counts and
180 sparse reference-Python scores, not all native scores. No coefficients,
thresholds or defaults were fitted or promoted.

### Evolutionary Simulations Retain Failures And Tree Sensitivity

In the fixed-length stress panel, high-sensitivity and phylogenetic OrthoHMM
had 70 and 64 admitted datasets out of 70. No OrthoFinder output passed its
native-completion/finite-graph gate; all 14 planned comparisons are unavailable,
not wins for OrthoHMM. The heterogeneous-length panel admitted 70 and 67
OrthoHMM datasets and 65 full OrthoFinder datasets. Full OrthoFinder had higher
paired mean F1 than both OrthoHMM modes in all seven conditions. Adjusted F1
intervals were below zero for all seven high-sensitivity contrasts and four
phylogenetic contrasts. Phylogenetic-minus-full-OrthoFinder effects were
-11.7433 points for divergence (five paired seeds) and -12.2842 for combined
divergence/turnover (eight). These are success-conditioned comparisons, not
failure-adjusted population estimates. Species-tree inference failures for
OrthoHMM and nonfinite graph weights for OrthoFinder remain explicit.
[Fixed stress results](SIMULATION_FIXED_NATIVE_INTERPRETATION_20260916.md),
[heterogeneous-length results](SIMULATION_VARIABLE_NATIVE_INTERPRETATION_20260916.md).

The complete tree-control analysis retained 560 arm outcomes: 537 scored and
23 failed. No generating-minus-inferred adjusted interval excluded zero. The
stronger perturbation reduced F1 and recall with adjusted intervals below zero
in four OrthoHMM conditions and two OrthoFinder conditions: 12 endpoints;
no precision interval excluded zero. Retained upstream artifacts agreed in
400 comparisons, differed in two and were unavailable in 18. The differing
OrthoFinder cases prevent a strict tree-only causal attribution there. These
deterministic perturbations and ten planned seeds do not establish robustness
to arbitrary trees or empirical posterior uncertainty.
[Complete tree results and figure](SIMULATION_TREE_ROBUSTNESS_RESULTS_20260917.md).

### Gene-Tree And Event-History Controls Localize Errors

Across 10,125 fixed candidates, 1,377 received eligible gene-tree controls;
8,748 bypasses remained unchanged. Among eligible controls, 300 differed in
unrooted topology and 30 additional candidates only in root. Generating-root
mean F1 changed by +0.136 to +0.734 percentage points across the seven
conditions. These descriptive effects are not significance tests or evidence
of independent improvement. In divergence and divergence/turnover, respectively
12,736 of 12,766 and 13,072 of 13,104 residual generating-root false negatives
(99.765% and 99.756%) already cross candidate boundaries. Changing gene trees
alone cannot recover them.
[All-condition oracle results](SIMULATION_GENE_TREE_ORACLE_RESULTS_20261004.md).

Of those cross-candidate losses, 11,442 and 11,554 (89.840% and 88.387%) lie
in different saved graph components. The remaining 1,294 and 1,518 are connected
but separated by candidates, including 117 and 128 direct graph edges.
Almost all component-separated losses lack significant initial hits in either
orientation. This identifies retained stage boundaries, not a single cause:
missing hits cannot isolate the prefilter, and paths may cross other families.
[Complete upstream trace](SIMULATION_UPSTREAM_TRACE_REVIEW_20261004.md).

The complete residual cohort contains 20 candidates: eight eligible controls
and twelve bypasses. All 800 cross-species pair rows reproduce 380 TP, 292 TN,
62 FN and 66 FP. All 62 false negatives receive high-confidence speciation
calls before unsupported satellite sources are detached; their root groups
agree, but the post-constraint filter removes them. Of the false positives,
45 are bypass calls at true duplication ancestors and 21 receive speciation
calls without retained species overlap. All 66 have overlap evidence in their
full parent generating family that is absent from their retained candidate.
Thus a single retained copy per species is not sufficient to establish
orthology. These selected-cohort counts do not replace whole-panel accuracy.
[Event-history residual results](SIMULATION_ORACLE_RESIDUAL_REVIEW_20261004.md).

The trace supports candidate-retention and membership-filter explanations for
these fixed simulation cases. It does not establish that disabling constraints,
merging candidates indiscriminately, or changing a search threshold improves
independent real-data inference. All negative and neutral findings, bypasses
and upstream losses remain. QfO/OrthoBench retain their primary role; these
privileged-information simulations are mechanistic support, not substitutes
for curated-reference validation or unfinished biological error strata.

### Local Parameter Sensitivity Did Not Establish Improvement

The [complete seven-arm panel](QFO_COMPLETE_PARAMETER_UNCERTAINTY_20261001.md)
includes separately admitted private-runtime high-CPM recovery. Its SwissTrees
F1 was 83.253016% versus 83.351322% for the frozen control, a difference of
-0.098306 percentage points with adjusted interval [-1.050447, 0.747525].
All 18 adjusted intervals include zero; no statistically supported parameter
improvement is established. Identical observed normalization-arm statistics
do not imply identical whole-proteome predictions or equivalence on unseen
families. No defaults were changed. The [all-arm figure](qfo_parameter_complete_export_20261001/qfo_parameter_neighborhood.pdf)
retains every contrast and the full adjustment denominator.

Historical SIGSEGV, allocator and admission failures remain failed. Recovery
followed a validated content-equivalent private control and has
[separate inference, conversion and score evidence](QFO_PRIVATE_CPM_SCORE_RESULT_22394.md);
it neither repairs those earlier attempts nor proves their cause or memory
safety. This panel is not an OrthoFinder superiority test or a controlled
runtime comparison.

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

### Shared-Host Resource Panel

The snapshot retains 27 of 27 reviewed attempts, 25 with measured native resources and 24 eligible observations. 6 method/size cells have three eligible repeats.

These are shared-host matched-resource observations with 32 native CPUs and a 128-GiB limit per run, not isolated comparative timing. Contention distortion is unknown and potentially method dependent. No background overhead is subtracted.

Tables report median [minimum, maximum] only for cells with three eligible repeats. An unavailable summary is not zero or a failed native inference. Ranges are not confidence intervals. Eligibility counts retain measured exclusions and pre-native aborts rather than selecting the fastest attempts.

#### Native command wall (seconds)

| Method | Proteomes | Eligible/planned | Median [minimum, maximum] |
| --- | --- | --- | --- |
| OrthoHMM high sensitivity | 4 | 2/3 | Unavailable |
| OrthoHMM high sensitivity | 8 | 3/3 | 1238.9352 [1221.4444, 1303.6330] |
| OrthoHMM high sensitivity | 12 | 3/3 | 2720.0774 [2621.4976, 3027.0605] |
| OrthoHMM inferred phylogeny | 4 | 2/3 | Unavailable |
| OrthoHMM inferred phylogeny | 8 | 3/3 | 1999.5508 [1888.8847, 2086.9473] |
| OrthoHMM inferred phylogeny | 12 | 3/3 | 4392.8716 [4042.4262, 4809.2924] |
| OrthoFinder 3.1.5 full | 4 | 3/3 | 441.2890 [420.7400, 474.2999] |
| OrthoFinder 3.1.5 full | 8 | 3/3 | 1139.9398 [1120.2548, 1179.9703] |
| OrthoFinder 3.1.5 full | 12 | 2/3 | Unavailable |

#### Task-subtree CPU bracket (CPU-seconds)

| Method | Proteomes | Eligible/planned | Median [minimum, maximum] |
| --- | --- | --- | --- |
| OrthoHMM high sensitivity | 4 | 2/3 | Unavailable |
| OrthoHMM high sensitivity | 8 | 3/3 | 34395.7481 [33878.9946, 36610.5696] |
| OrthoHMM high sensitivity | 12 | 3/3 | 78335.8306 [75126.8588, 86415.3192] |
| OrthoHMM inferred phylogeny | 4 | 2/3 | Unavailable |
| OrthoHMM inferred phylogeny | 8 | 3/3 | 56008.3471 [52911.7819, 58240.7234] |
| OrthoHMM inferred phylogeny | 12 | 3/3 | 125264.7375 [116000.9964, 136293.3322] |
| OrthoFinder 3.1.5 full | 4 | 3/3 | 4602.7585 [4455.3737, 5008.8970] |
| OrthoFinder 3.1.5 full | 8 | 3/3 | 17734.4088 [17232.2843, 18787.3118] |
| OrthoFinder 3.1.5 full | 12 | 2/3 | Unavailable |

#### Native-step lifetime peak (GiB)

| Method | Proteomes | Eligible/planned | Median [minimum, maximum] |
| --- | --- | --- | --- |
| OrthoHMM high sensitivity | 4 | 2/3 | Unavailable |
| OrthoHMM high sensitivity | 8 | 3/3 | 7.724 [7.715, 7.932] |
| OrthoHMM high sensitivity | 12 | 3/3 | 10.555 [10.533, 10.730] |
| OrthoHMM inferred phylogeny | 4 | 2/3 | Unavailable |
| OrthoHMM inferred phylogeny | 8 | 3/3 | 7.681 [7.677, 7.702] |
| OrthoHMM inferred phylogeny | 12 | 3/3 | 10.511 [10.503, 10.603] |
| OrthoFinder 3.1.5 full | 4 | 3/3 | 6.775 [6.769, 6.928] |
| OrthoFinder 3.1.5 full | 8 | 3/3 | 10.226 [10.215, 10.233] |
| OrthoFinder 3.1.5 full | 12 | 2/3 | Unavailable |

Excluded attempt indices: 0, 17, 20. Pre-native abort indices: 17, 20; these have no native resource endpoints, not zero measurements. Unreviewed indices: none.

Across measured attempts, maximum observed foreign CPU demand ranges from 42.1424 to 72.8708 core equivalents. These are process-interval observations, not reservations or estimates of causal slowdown.

CPU includes the native task-subtree wrapper bracket. Peak memory includes the native-step launcher and is not pure algorithm RSS. Preparation, conversion and scoring are outside the native command timer. This single nested taxon series co-varies proteome count and taxon composition, not taxon-invariant scaling. Native output validation does not provide new prediction-accuracy evidence. Historical DGX and earlier shared-host timings are not pooled.

The first high-sensitivity attempt completed native inference but failed the
original process-monitoring cadence criterion. It remains a measured exclusion;
the prospectively tested scheduling repair does not convert it to a pass.
Two further attempts aborted before native inference: an OrthoFinder/twelve-
proteome repeat encountered a preparation hash race, and a phylogenetic
OrthoHMM/four-proteome repeat exceeded the environmental-response deadline.
Both retain null native endpoints and separately reviewed infrastructure
repairs. None of these failures was retried to fill a cell or obtain faster
timing. All 27 planned identities are terminal-reviewed, including exclusions;
this is not 27 eligible measurements. The three incomplete cells have two
eligible repeats each, not completed three-repeat summaries.

[Machine-readable final table](threadripper_shared_panel_snapshot_20261004_v27/panel.json),
[all attempts](threadripper_shared_panel_snapshot_20261004_v27/attempts.tsv),
[complete-cell summaries](threadripper_shared_panel_snapshot_20261004_v27/cells.tsv),
[final resource figure](threadripper_shared_resource_figure_20261004_v26/shared_threadripper_resources.pdf),
[final terminal review](THREADRIPPER_SHARED_ATTEMPT_22424.md),
[generated section and provenance](threadripper_shared_resource_section_20261004_v27/manifest.json).

### Retained Development Exposure And Candidate Reproducibility

The family-evidence inventory contains 134 OrthoBench score blocks across
84 files and 77 SwissTrees blocks across 16 files: 7,770 and 1,386 family-block
associations respectively. All 88 canonical families have explicit score
evidence. Each of the 35 originally validation-labelled RefOGs appears in
18 declared validation blocks; each original development RefOG appears in
28 development blocks, and all 70 appear in 21 all-partition blocks. The
historical 35/35 split remains recorded chronology, not untouched publication
validation. Distinct count vectors number 84 and 36; neither different nor
equal vectors establish independence or duplicate execution. The inventory
also retains 34 subset-ordinal candidate-diagnostic containers and two
synthetic scored-label containers without misassigning biological identities.
[Verified inventory](DEVELOPMENT_FAMILY_INVENTORY_RESULT_20261004.md).

Each original/native candidate trace contains 8,440 accepted merges. The
three native repeats share 8,428, 8,439 and 8,426 semantic records with the
original. Common cluster indices and hit counts agree; the maximum common
support difference is 2.1316282072803006e-14. Every common anchor with changed
selected source memberships reaches the four-attachment cap in both traces.
In the frozen 19-gene fixture, baseline, order-reversed and one-unit-perturbed
cases leave satellite IDs 8, 0 and 7 unattached, respectively: one unattached
satellite and eight accepted merges in every case. These IDs are not numbers
of missing satellites. This demonstrates that capped selection
and tied/near-tied evidence can change assignments without changing merge
counts. It does not establish historical relabeling, score-bit provenance,
gene-tree correctness or an accuracy advantage.
[Complete diagnostic](CANDIDATE_TRACE_VARIATION_RESULT_20261004.md).

### Historical Costs Remain Scope-Specific

A retained native-cost linkage checks six completed twelve-proteome
Threadripper repeats against the original OrthoBench inputs and prescribed
settings. All three high-sensitivity final partitions equal the original
`p1_c0_r0` partition. The first phylogenetic final partition matches the
original `p1_c1_r1`; the other two have changed groups containing 149 and
172 genes. All three candidate partitions differ, involving 70, 160 and
253 genes in affected groups. No changed group intersects a reference family.
Affected-group gene counts are not counts of genes moved, and structural
agreement is not new accuracy scoring or equivalence of all intermediate
artifacts. Native wall medians are 45.335 and 73.215 minutes; these reuse the
shared-host panel rather than provide new measurements or costs of the original
cached executions. The patched native runtime differs from that deployment.
[Native cost/partition linkage](FACTORIAL_NATIVE_RESOURCE_LINKAGE_RESULT_20261004.md).

The selected corrected QfO high-sensitivity and phylogenetic rows bind to
`p1_c0_r0` group-derived pairs and `p1_c1_r1` native pairs respectively.
Their retained resource associations describe one shared cached-replay worker
(2,987.81 seconds), two preparation arms (3.035759820 and 77.641008615 seconds)
and native reconciliation (6,822.365357 seconds). The repeated replay value
is one observation, not two repeats; process RSS, sampled tree RSS and
lifetime cgroup peak retain different meanings. Neither stage sums nor these
native configuration observations fill the sixteen unavailable original
factorial full-pipeline costs. The selected rows' complete cached-execution
costs remain unknown.
[Exact QfO stage chains](QFO_ORTHOHMM_STAGE_LINKAGE_RESULT_20261004.md),
[factorial stage costs](FACTORIAL_RESOURCE_RESULT_20261004.md).

The additive all-tool metadata register retains 24 method/dataset rows and
all 72 scientific score positions, with explicit conversion semantics and
historical resource supplements. Selected QfO stage associations are added
without replacing prior intervals; table entries are not independent timed
experiments. Input consumption, transitive runtime attestation and missing
whole-pipeline costs remain separately identified rather than inferred.
[All-tool metadata integration](ALL_BENCHMARK_METADATA_INTEGRATION_RESULT_20261004.md).

## Discussion And Limitations

The supported contribution is an HMM-centered alternative with measurable
component effects and explicit accuracy trade-offs. Initial HMM search has
bounded support in matched-recall simulations; additional profile refinement
has not demonstrated a general benefit. Neither simulation evidence nor an
OrthoBench point advantage establishes superiority over full OrthoFinder.

The approximate significance formula is not universally calibrated under
the tested synthetic conditions. The null experiment does not measure real-data
orthology false-positive rates or causally explain any benchmark difference;
independent calibration would be needed before claiming statistical significance
for arbitrary compositions, search spaces, profiles or band settings.

Original TreeFam-A family mappings and complete source trees remain unavailable.
Public QfO container searches now cover recognized benchmark build contexts
for [all 40 listed tag names](TREEFAM_REMAINING_CONTEXT_RESULT_20260930.md)
at retained digests, without finding original NHX/mapping filename candidates.
Other layers and archive-embedded or differently named contents remain outside
that search; no family labels are inferred from its negative result.
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
The [requested-score audit](QFO_FAS_SAMPLE_ATTRITION_20260928.md) found 1-49
missing new scores per method except OrthoMCL, which lacked 1,252 of 9,000.
Assuming omitted scores lie in [0,1], its intended-sample mean is bounded by
0.724872-0.736975; this is not a confidence interval. A
[native-complexity audit](FAS_SAVED_COMPLEXITY_EXPOSURE_20260928.md) found no
above-limit protein in saved new-score pairs, but found them in precomputed
pairs. This agrees with the demonstrated cutoff mechanism without identifying
historical omitted pairs or establishing comparison bias.
The [completed eligible-population audit](QFO_FAS_POPULATION_COMPLETED_22383.md)
matches all eight methods' native logged counts, reusing six fixed recounts
and freshly scanning two. Pairs absent from the precomputed lookup constitute
5.13% of satellite_v2's eligible set, 16.86% of full OrthoFinder's and 82.05%
of its sequence-only checkpoint's. Assuming their scores lie in [0,1] gives
conservative full-mean bounds of 0.733966-0.785282 for satellite_v2 and
0.605946-0.774511 for full OrthoFinder. These overlapping bounds neither
establish a full-eligible-set advantage nor replace native sample scores or
confidence intervals. Some lookup-missing pairs were newly scored in native
samples; those values are not used to tighten these bounds.
An [eight-method stratum-weight decomposition](QFO_FAS_STRATUM_WEIGHT_RESULT_20261001.md)
separates count rounding from omission-induced changes in the saved mixture.
Holding the observed stratum means fixed, native-minus-reweighted-saved-strata
differences are +0.001318947 FAS units for OrthoMCL and +0.000030716 for
phylogenetic OrthoHMM. These are numerical diagnostics, not corrected benchmark
scores or estimates of population bias. Omissions can also change the stratum
means; the unknown omitted pair identities, unseeded selection and dependence
remain unresolved. No FAS confidence interval or new method ranking follows.
The completed simulations and tree perturbations do not cover arbitrary
evolutionary conditions, and novel-taxon YGOB testing retains homolog-family
overlap with development data.
Generating-tree, upstream and residual diagnostics localize errors in the
retained development-exposed panel but do not establish empirical gene-tree
accuracy, ancestral-copy root-HOG truth or independent generalization.
The development inventory establishes retained evaluation exposure, not a
complete causal tuning history, all historical experiments or family-disjoint
validation. Its frozen local input set must be supplied to reproduce that
expanded inventory. Capped tie-sensitive candidate selection is a demonstrated
reproducibility limitation, but does not justify benchmark-chosen score rounding
or attributing changes to competing workloads. A method change would require
a new freeze and independent confirmation. Candidate retention can remove species-overlap evidence, and strict satellite detachment
can remove otherwise correct pair calls. The demonstrated cases motivate
testable hypotheses, not unvalidated new defaults.

Recovered OrthoMCL results retain sequence-specific BLAST failures. The
separately recovered high-CPM result does not replace failed historical runs;
no missing score was imputed and no failure was recategorized as a success.
Shared-host runs and the historical DGX panel are descriptive resource records.
The replacement matched-resource Threadripper panel has completed under the
shared-host amendment, not as an isolation-controlled experiment. Competing
workloads may affect tools differently, and identical resource limits do not
identify isolated speed or causal speedups. The final figure retains excluded
measurements and missing eligible repeats; it cannot support an isolated
efficiency ranking. Six cells have three eligible repeats, while three cells
retain unavailable summaries. This single nested taxon series does not isolate
taxon-count scaling from proteome composition. Whole-study executable/versioned
archival packaging remains open, alongside redistribution limitations and
journal-specific formatting. Completed resource reporting does not resolve the
scientific uncertainty and validation gaps above.

## Reproducibility And Availability

The extended manuscript links frozen protocols, source revisions, seeds,
input/output manifests, conversion audits, executable workflows and generated
figures. Installed inference and independent readback have been validated on
full OrthoBench, but that is not proof of all-method cross-host portability.
Historical locks are retained as provenance; patched replacement environments
must be distinguished from the binaries used for original scores. The
[QfO parameter numerical component](QFO_PARAMETER_NUMERICAL_COMPONENT_RESULT_20261001.md)
passed fresh local archive extraction, copied byte verification and reproduction
of all 18 endpoints within 1e-12 using unchanged arithmetic. A Python-event
guard rejected its original-checkout canary and observed no later original-path
events or project-module imports. This is not OS containment, a new independent
statistical implementation, raw-reference recount, native inference restoration,
cross-host validation or redistribution clearance.

A separate [descriptive-table component](SWISS_DESCRIPTIVE_COMPONENT_20261002.md)
restored under standard-library Python 3.12 and reproduced all 208 rows and
984 score/difference cells across four SwissTrees tables and eight methods.
This verifies arithmetic from retained sufficient statistics and annotations,
not raw-source admission, new benchmark scoring or bootstrap uncertainty.
An [input-only private archive workflow](SWISS_RAW_ARCHIVE_RESTORATION_20261002.md)
then restored the complete duplication and fragment dependency panels, retaining
all 11 and 1,774 original record occurrences respectively. Independently pinned
archive and binding digests preserve source identities at new locations.
All 27 affected regression/options cases passed using the restored inputs,
including four raw-export regressions with unchanged assertions. The fresh
children blocked original-checkout reads and subprocesses; these Python-event
guards are not OS containment. Native annotation extraction/admission was not
rerun, and private archives remain unuploaded with redistribution uncleared.
Neither component establishes complete executable study restoration or
all-method cross-host portability.

The [null-score observations](frozen_null_score_observations_20261002.json.gz)
retain all 180,000 synthetic score evaluations and their source/runtime pins.
Its public numerical readback checks all 90 tail endpoints without loading
the native kernel; this is not a complete inference runtime archive.

The [standalone YGOB arithmetic replay](YGOB_ARITHMETIC_REPLAY_RESULT_20261002.md)
uses identifier-free retained counts to reproduce all 12 point metrics and
six nominal/adjusted interval endpoints after local archive restoration.
This is not native re-admission, proof of pillar exchangeability or independent
family validation; original transfer scores and overlap limitations remain.

The [standalone simulation replay](../SIMULATION_ARITHMETIC_REPLAY.md) reproduces
retained-count arithmetic for the two length panels without historical file
reads. It preserves fixed-length unavailable comparisons and variable-length
conditional effects; native admission and tree-control inference are not rerun.

The [final resource-reporting archive](FINAL_RESOURCE_REPORTING_COMPONENT_20261004.md)
restores outside the repository and reproduces all three table files, PNG
pixels and resource-section text exactly from its copied isolated reader.
It includes all attempts and failures without opening original evidence paths.
This is reporting portability, not raw-accounting or native reproduction,
isolated timing, complete study restoration or redistribution clearance.
The [source/review export integration](PUBLICATION_EXPORT_INTEGRATION_20261004.md)
adds explicit result-helper and final-manuscript selection; fixture tests do
not establish an actual full-study archive.

The
[progress ledger](PUBLICATION_PROGRESS.md) records completed work and unmet
requirements. No submission-ready release or archival DOI is claimed.
Dated main-text HTML/PDF and review archives preserve the exact inputs of
their own snapshots; older versions are not rendered copies of revised
Markdown. A refreshed review is not a complete executable study release.
This revision incorporates the completed shared-host resource panel while
also integrating the later simulation mechanism diagnostics and preserving
scientific limitations. The three hash-pinned checked-summary reports support
a standalone reporting replay using standard-library Python; it recomputes
the finite-panel means, stage-count identities and residual classifications,
not native search/inference or original XML admission. Its separate execution
and package receipts govern what has actually been restored and reproduced.
The second revision's 34-page review/rc3 archive and the third revision's
35-page review/rc4 archive remain historical. The latter includes the exposure,
reproducibility and cost paragraphs, but neither archive includes this
6 October native-evidence revision. Fresh native scientific tables, bounded
interval bindings, the reviewed QfO figure and functional-pair diagnostic
retain their own source/admission identities and failure records.
The new source revision requires its own render/review and eventual versioned
bundle integration; those steps have not been executed here. Existing native
and statistical replay receipts remain valid only within their recorded scopes.
No new endpoint selection, score, method default or independent-validation
claim accompanies this revision.

## References

The bibliography attributes the methods and resources cited above. Full
execution evidence and additional dependency references remain in the extended
manuscript; this selected list is not complete software attribution.
