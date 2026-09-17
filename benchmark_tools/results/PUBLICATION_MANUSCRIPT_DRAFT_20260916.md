# OrthoHMM: HMM-Centered Orthogroup Inference And Phylogenetic Refinement

Working manuscript, 16 September 2026. Not submission-ready. Sections below
distinguish completed development-exposed analyses from prospective work.
Internal evidence links are supplied for audit; a verified literature
bibliography and journal-specific formatting remain to be added.

## Study Objective

We evaluate an HMM-centered approach to orthogroup inference and its
integration with phylogenetic refinement. The principal questions are whether
profile-based expansion contributes measurable recovery beyond the initial
search, whether broader candidate groups improve reconciliation, and how
accuracy and computational costs compare with established tools. These
questions are not equivalent to asking whether a single aggregate score is
maximized. We distinguish homolog-group recovery, resolved ortholog-pair
prediction, and conserved-family recovery throughout the study.

## Methods

### Configurations And Comparator Outputs

The retained OrthoHMM configurations are high sensitivity and the
satellite_v2 phylogenetic pipeline. The prospective validation configuration
is pinned to source revision `7f3a9e4`, with BLOSUM62, E-value threshold
1e-4, Leiden CPM resolution 0.1, seed 4, and the default refinement profile.
The phylogenetic configuration uses satellite_v2 candidate expansion,
an internally inferred species tree, minimum-variance species-tree rooting,
the species-overlap root rule, and positive-paralogy pair inference.
Historical runs retain their actual source and configuration records; the
prospective pin must not be retroactively attributed to them.
See the [frozen validation protocol](YGOB_VALIDATION_PROTOCOL_20260916.md).

The [method diagram](figures_publication_method_20260916/publication_method.pdf)
distinguishes initial HMM search, cluster-profile expansion, high-sensitivity
orthogroups, candidate-family construction, gene/species-tree inference and
constraint-guided membership decisions. Root HOGs and inferred ortholog pairs
are separate output levels. This schematic describes the intended frozen
workflow, not evidence of execution correctness or an HMM accuracy advantage.

Comparators comprise OrthoFinder 3.1.5 full and its sequence-only MCL
checkpoint, OrthoMCL 1.4, SonicParanoid 2.0.9, ProteinOrtho 6.3.6, and
FastOMA 0.3.5. The checkpoint is a distinct output level, not another
phylogenetically finalized prediction. In the retained OrthoFinder full
run, finalized `Orthogroups.txt` reflects root hierarchical-group
postprocessing; its filename does not identify it as sequence-only.
FastOMA's supplied OrthoFinder tree is disclosed rather than treated as
independent tree inference. Native-output semantics and source records are
retained in the [comparison](PUBLICATION_COMPARISON_ORTHOMCL_COMPLETE_20260916.md)
and [OrthoBench protocol](ORTHOBENCH_UNCERTAINTY_PROTOCOL_20260916.md).

For QfO, OrthoHMM high sensitivity and the OrthoFinder checkpoint contribute
group-derived cross-species pairs; OrthoHMM satellite_v2 and full OrthoFinder
contribute phylogenetically inferred pairs. SonicParanoid contributes native
species-pair relations, ProteinOrtho its native post-clustering graph, and
FastOMA native pairs. OrthoMCL contributes cross-species pairs derived from
its final MCL groups. Pre-MCL graph edges are retained only as a diagnostic.
These output differences are part of the practical pipeline comparison and
must not be interpreted as a matched search-engine ablation.

### Development-Exposed Benchmarks

QfO and OrthoBench influenced development and are therefore not independent
confirmation sets. OrthoBench evaluates the complete retained 70-RefOG panel
with its low-certainty exclusions and audited weighted pair statistic.
Weighted counts divide each reference-family contribution by reference size
minus one; aggregate precision, recall, and F1 are computed from those
counts, not by averaging family F1 values.

We resampled RefOGs jointly across methods 20,000 times using NumPy PCG64
seed 20260916 and recomputed the statistic for each replicate. We report
paired percentile intervals and Bonferroni-adjusted intervals across two
OrthoHMM-versus-OrthoFinder contrasts and three metrics. The protocol was
recorded after historical aggregate results were known, before this
bootstrap; it is not preregistration of method selection. Family
exchangeability and residual cross-family dependence limit inference.
[Protocol](ORTHOBENCH_UNCERTAINTY_PROTOCOL_20260916.md).

For QfO we retain individual VGNC, SwissTrees, TreeFam-A, EC, GO, and FAS
endpoints, their native axes, and source provenance. The unweighted mean of
six project-selected summaries is secondary and is not an official QfO
score. Relation counts assessed by functional endpoints do not necessarily
equal total prediction coverage. Three Kingdoms is supplementary: its
BUSCO-reference pair statistic ignores false positives involving
non-reference genes and cannot establish proteome-wide orthology accuracy.
[Machine-readable comparison](publication_comparison_orthomcl_complete_20260916.json).

### Prospective Validation And Ablations

The frozen YGOB experiment retains 16 non-Saccharomyces species, 83,404
inference proteins, and 83,391 reference genes in 10,250 curated pillars.
Thirteen proteins from ambiguous duplicate-membership rows remain in
inference but are projected out of scoring. The endpoint is homolog-group
co-membership, including within-species pairs, not resolved orthology.
The primary contrast is satellite_v2 versus full OrthoFinder; high
sensitivity is secondary and the sequence checkpoint diagnostic. Its
20,000-replicate paired pillar bootstrap has seed 20260917 and six-metric
multiplicity handling. Native/source/input, command/conversion, independent
reference reconstruction and overlap checks passed before scoring. Direct
pair enumeration independently reproduced all four methods' TP/FP/FN counts.
[Frozen protocol](YGOB_VALIDATION_PROTOCOL_20260916.md).

On this frozen panel, satellite_v2 achieved F1 92.233654%, versus 92.318524%
for full OrthoFinder. The paired difference was -0.084870 percentage points
(nominal 95% interval [-0.487503, 0.312676]; six-endpoint Bonferroni interval
[-0.622528, 0.445222]). This does not establish superiority or equivalence.
Satellite_v2 had higher precision (+6.401035 points) and lower recall
(-6.914033 points), with both adjusted intervals excluding zero. High
sensitivity scored 82.038608% F1, below full OrthoFinder by 10.279916 points
(adjusted interval [-11.260924, -9.298064]). The diagnostic sequence-only
OrthoFinder checkpoint scored 85.572844%. All methods covered every reference
gene. No configurations were changed in response to these held-out outcomes.
[Frozen results](YGOB_FROZEN_RESULTS_20260916.md),
[machine-readable evidence](ygob_frozen_results_20260916.json).

![Frozen YGOB scores and paired differences](figures_ygob_frozen_20260916/ygob_validation.png)

Figure: observed co-membership scores and paired OrthoHMM-minus-OrthoFinder
differences. Thick intervals are nominal 95%; thin intervals apply the frozen
six-endpoint Bonferroni adjustment. The checkpoint is diagnostic. All methods
cover every scored reference gene; this does not establish family independence.

This experiment tests novel-taxon transfer, not family-disjoint validation.
The completed label-independent homology screen found qualifying development
input hits for 71,714 of 83,404 proteins and 6,952 of 10,250 pillars.
Absence of a qualifying hit does not establish absence of homology.
Shared reference-resource limitations are documented separately.
[Screen](ygob_homology_screen_20260916.json),
[resource audit](YGOB_REFERENCE_RESOURCE_AUDIT_20260916.md).

Prospective component experiments cross profile-HMM expansion, candidate
expansion, and reconciliation in eight cells. Expansion-plus-reconciliation
cells retain their own membership constraints; independently inferred trees
are distinguished from supplied-tree diagnostics. A profile-expansion-off
arm still uses the HMM-centered initial search. A matched sequence-search
control is required before attributing an overall advantage to HMMs.
The OrthoBench factorial is complete; its QfO counterpart and additional
matched-search and membership-filter controls remain unfinished.
[Ablation protocol](PUBLICATION_ABLATION_PROTOCOL_20260916.md).

### Evolutionary Simulations And Runtime Admission

Two separate panels each contain ten independent simulation seeds and seven
conditions: baseline, increased divergence, increased duplication/loss,
combined divergence and turnover, missing proteins, uneven taxon sampling,
and a taxon-count control. Native Zombi histories provide event-based
cross-species orthology truth; shared ancestral-family membership alone is
not treated as resolved orthology. Four configurations are evaluated on each
of the 70 datasets. The first panel uses fixed 300-residue sequences; a
prospectively specified second panel assigns each ancestral family a
deterministic seed-specific length between 100 and 500 residues. Descendants
retain that length, so this panel does not simulate within-family indels.
[Fixed protocol](PUBLICATION_SIMULATION_PROTOCOL_20260916.md),
[variable-length protocol](PUBLICATION_VARIABLE_LENGTH_PROTOCOL_20260916.md).

Initial OrthoHMM runs lacked a required native alignment library; profile
construction exceptions were silently caught and no profiles were built.
Those runs are invalid runtime diagnostics, not estimates of the intended
method. Corrected runs retain the frozen algorithm and verify compiled
libraries and an actual profile-construction probe before and after execution.
Valid original comparator outputs are reused with separate provenance.
Admission requires native completion and output validation, including finite
OrthoFinder graph weights; a zero process exit code alone is insufficient.
[Runtime correction](SIMULATION_RUNTIME_CORRECTION_PROTOCOL_20260916.md).

Within each panel, we bootstrap paired successful seeds 20,000 times and
recompute the mean seed-level statistic. Fourteen primary F1 contrasts use
Bonferroni-adjusted intervals; precision and recall are exploratory. Failed
runs are reported separately, not assigned zero accuracy. Complete-case
contrasts condition on joint success, and available-case table means cannot
be subtracted when they summarize different seeds. The two panels are not
pooled. [Fixed results](SIMULATION_FIXED_NATIVE_RESULTS_20260916.md),
[variable results](SIMULATION_VARIABLE_NATIVE_RESULTS_20260916.md).

## Results

### OrthoBench Shows A Precision-Recall Tradeoff

OrthoHMM satellite_v2 achieved weighted F1 74.106074%, compared with
72.736480% for full OrthoFinder and 70.358998% for high sensitivity.
Satellite_v2 precision was 81.770454% and recall 67.755336%; full
OrthoFinder precision was 66.065103% and recall 80.906577%.
The satellite_v2 F1 difference was +1.370 percentage points, with nominal
95% interval [-4.504, 7.917] and adjusted interval [-6.290, 10.627].
These intervals do not establish an F1 advantage. Its precision difference
was +15.705 points and recall difference -13.151 points; corresponding
adjusted intervals excluded zero in opposite directions. These are
development-exposed comparisons, not selection-adjusted confirmation.
[Uncertainty results](ORTHOBENCH_UNCERTAINTY_20260916.md).

Across individual RefOGs, satellite_v2 had 23 F1 wins, 10 ties, and 37
losses relative to full OrthoFinder. High sensitivity had 20 wins, 11 ties,
and 39 losses. These descriptive counts need not rank methods identically
to the weighted aggregate statistic and do not substitute for it.
[Per-family analysis](orthobench_paired_uncertainty_20260916.json).

### QfO Results Vary Across Endpoints

The secondary QfO means were 0.782071 for full OrthoFinder, 0.748243 for
satellite_v2, and 0.682548 for high sensitivity. Full OrthoFinder had higher
VGNC, SwissTrees, and TreeFam-A F summaries than either OrthoHMM mode.
Satellite_v2 had higher EC, GO, and FAS similarity summaries than full
OrthoFinder. Because these functional summaries assess different predicted
relation sets, higher similarity alone does not establish higher overall
orthology accuracy. The six native endpoint panels retain their assessed
relation counts or precision-recall axes.
[Comparison](PUBLICATION_COMPARISON_ORTHOMCL_COMPLETE_20260916.md),
[endpoint figure](figures_accuracy_orthomcl_complete_20260916/qfo_endpoints.png).

OrthoMCL final-group scoring completed successfully, with all six endpoint
tasks and consolidation exiting zero. Its secondary mean was 0.705177.
The earlier pre-clustering diagnostic mean, 0.724414, is not the final-group
assessment and is excluded from the principal comparison. The run contained
53 sequence-specific BLAST failures among 976,504 proteins; all 53 were
absent from final groups. Direct-reference exposure was absent for VGNC,
SwissTrees, TreeFam-A, and EC, but present for four experimentally annotated
GO proteins; 46 failed proteins had nonempty FAS features. These observations
do not bound indirect clustering effects or establish a negligible score
impact. The unmodified baseline is retained with failures disclosed.
[Failure-impact audit](ORTHOMCL_FAILURE_IMPACT_20260916.md).

### Conserved-Family Recovery Remains A Weakness

Three Kingdoms BUSCO-reference F1 was 0.872133 for satellite_v2 and
0.826309 for high sensitivity, compared with 0.988582 for full OrthoFinder,
0.989451 for its sequence checkpoint, and 0.990794 for SonicParanoid.
The checkpoint's small advantage over full OrthoFinder on this restricted
target does not establish that phylogenetic inference is generally harmful.
Mechanistic explanations require error tracing and controlled ablations.
[Comparison](PUBLICATION_COMPARISON_ORTHOMCL_COMPLETE_20260916.md).

### Historical Components Provide Preliminary Evidence

Verified historical OrthoBench partitions gave F1 66.279051% after
multipass grouping, 69.763388% with cluster refinement, 66.825390% with
profile expansion alone, and 70.358998% with both. The descriptive increment
from profile expansion after refinement was +0.595610 percentage points.
The final replay partition was byte-identical to the corrected historical
fresh output. These are historical component observations without paired
intervals, not a complete current-source factorial or HMM-free control.
[Historical audit](HISTORICAL_PROFILE_ABLATION_20260916.md).

### Corrected Simulations Do Not Establish An OrthoHMM Advantage

In the fixed-length stress panel, high sensitivity completed all 70 datasets
and satellite_v2 completed 64. Six satellite_v2 failures persisted because
species-tree inference lacked connected single-copy taxon coverage. No full
OrthoFinder or parent-gated checkpoint output passed native admission, so
this panel provides no admitted paired comparison with OrthoFinder. Native
normalization diagnostics reproduced degenerate length fitting and nonfinite
scores; competitor failure is not evidence of an OrthoHMM accuracy advantage.
[Fixed interpretation](SIMULATION_FIXED_NATIVE_INTERPRETATION_20260916.md),
[native failure audit](SIMULATION_NATIVE_FAILURE_AUDIT_20260916.md).

In the variable-length panel, high sensitivity, satellite_v2, and full
OrthoFinder passed admission on 70, 67, and 65 datasets, respectively. Full
OrthoFinder led every paired condition-mean F1 comparison against either
OrthoHMM mode. Adjusted intervals excluded zero below OrthoFinder for all
seven high-sensitivity contrasts and four satellite_v2 contrasts; satellite_v2
turnover, missing-data, and uneven-sampling intervals included zero. Baseline
satellite_v2 and full OrthoFinder mean F1 were 99.46% and 99.94%. Satellite_v2
paired deficits were 11.74 percentage points under divergence (five paired
seeds) and 12.28 points under combined divergence and turnover (eight seeds).
Lower recall was the principal observed deficit, but its causal processing
stage has not yet been established.
[Variable interpretation](SIMULATION_VARIABLE_NATIVE_INTERPRETATION_20260916.md).

Three variable-panel satellite_v2 failures again involved species-tree
coverage. Five divergent OrthoFinder runs had nonfinite graph weights.
Reconstruction of native normalization across all ten divergent seeds
reproduced seven nonfinite within-species matrices in exactly those five
failed seeds. Each affected fit had two non-self hits at one length product;
the rank-deficient fit produced extreme intercepts and scale-factor overflow.
The proteomes themselves had heterogeneous lengths. This local failure does
not imply that every rank-deficient fit fails, and no competitor repair or
post-hoc exclusion of additional finite runs was applied.
[Normalization audit](ORTHOFINDER_VARIABLE_NORMALIZATION_AUDIT_20260916.md).

All corrected OrthoHMM runs constructed profiles, but neither simulation
panel added profile-expansion edges. These panels therefore do not demonstrate
an accuracy contribution from multi-sequence profile expansion, although the
initial HMM-based search remains present. They are not HMM-free controls.
Synthetic sequence evolution, limited seed counts, conditional comparisons,
and absent realistic domain architecture constrain extrapolation to proteomes.

The [variable-length figure](figures_simulation_variable_native_v2_20260916/simulation_evidence.png)
shows admitted-run counts beside paired F1 effects and nominal/adjusted
intervals. The [fixed-length figure](figures_simulation_fixed_native_v2_20260916/simulation_evidence.png)
explicitly shows the absence of admitted comparator contrasts. Both are
generated from the corrected machine-readable results without imputing failures;
PDF/SVG versions and source-hash manifests accompany the PNGs.

### Controlled Components Show A Candidate-Precision Tradeoff

The completed eight-cell OrthoBench factorial reproduced full-pipeline F1
74.106074%. Reconciliation improved F1 by 2.942-6.899 percentage points across
four matched settings; all four multiplicity-adjusted intervals excluded zero
above it. Broader candidates increased recall and reduced precision, with
adjusted intervals excluding zero in all four settings. Their observed F1
effect changed from about -3.1 points without reconciliation to +0.697/+0.795
with reconciliation, but every adjusted candidate-expansion F1 interval
included zero. Profile expansion added 0.568-0.704 observed F1 points;
all adjusted profile-expansion intervals included zero. These contrasts do
not remove the initial HMM search and do not establish a total HMM advantage.

All cells passed output/provenance gates and official-score crosschecks.
Four batch failures caused by postflight working-directory-dependent package
discovery were retained and separately recovered; native inference succeeded.
Shared-node cached reconciliation costs are not end-to-end timings. The
factorial remains development-exposed and supplies no new direct OrthoFinder
comparison. [Results and interpretation](ORTHOBENCH_FACTORIAL_INTERPRETATION_20260916.md).

![OrthoBench factorial scores and paired conditional effects](figures_orthobench_factorial_20260916/orthobench_factorial.png)

Figure: all eight cell scores and 36 conditional-effect endpoints. Thick
intervals are nominal 95% paired RefOG bootstrap intervals; thin intervals
retain Bonferroni adjustment across all 36 endpoints. Each effect holds the
other two factors fixed. Profile-off retains initial HMM search; this figure
does not compare an HMM-free method or establish independent superiority.

An exploratory unconstrained replay used the same expanded candidate families
and reconciliation rules but omitted membership constraints. Native output
validation and official-score crosschecks passed. Its F1 was72.614510%,
precision76.092405%, and recall69.440641%, compared with74.106074%,81.770454%,
and67.755336% under constraints. Unconstrained-minus-constrained F1 was
-1.492 percentage points (three-endpoint Bonferroni interval[-5.873,0.931]);
precision decreased5.678 points ([-14.307,-0.794]) and recall increased1.685
points ([0.013,3.804]). Family F1 improved in5, tied in55, and declined in10
RefOGs. These development-exposed results describe a precision-recall tradeoff,
not a demonstrated F1 improvement or statistically established F1 loss.
They do not justify removing constraints as a new default.
[Diagnostic results](ORTHOBENCH_UNCONSTRAINED_RESULTS_20260916.md).

### Sequence-Search Replacement Control

With profile expansion, candidate expansion and reconciliation disabled,
the frozen HMM-search baseline scored69.763388% F1,78.868592% precision and
62.542944% recall. DIAMOND all-hit replacement scored65.762183/55.036231/
81.680884%; the post-search top100 diagnostic scored66.484803/55.897861/
82.019046%. Six-endpoint adjusted F1 intervals for control-minus-HMM were
[-14.834,6.688] and[-14.424,7.844] percentage points. Both include zero;
these results do not establish an HMM-specific F1 advantage. Precision was
lower and recall higher for both controls, with adjusted intervals excluding
zero. Each control improved family F1 in39, tied in10, and declined in21
RefOGs; the aggregate is not a mean of those family scores.

Hit sets were not sensitivity-matched. HMM/all-hit/top100 nonself counts were
18,235,373/99,848,167/48,740,683. Queries lacking cross-species hits numbered
33,351/24,526/24,526. DIAMOND all-hit recovery overlapped65.12% of HMM hits,
while only11.89% of its nonself hits were shared with HMM. This asymmetry is
not ground-truth sensitivity. The historical HMM cache omits self hits,
whereas DIAMOND retains250,980; raw no-hit query counts are not commensurate.
All three final partitions cover251,378 input genes, including singletons.
The same normalization formula and E-value cutoff do not equate calibration,
candidate filtering or computational effort. [Full control results](OB_SEQUENCE_SEARCH_RESULTS_20260916.md).

[Sequence-search control figure](figures_ob_sequence_search_20260916/sequence_search_control.pdf)
shows all observed scores and all six paired effects on a common percentage-point
axis. Thin intervals retain the six-endpoint correction; thick intervals are
nominal95% intervals. Search sensitivity/calibration are not matched, and this
control disables subsequent profile expansion and phylogenetic inference.

### Supplied-Tree Perturbations Had Small Observed Aggregate Effects

Six prespecified rooted topology perturbations were evaluated against an
unchanged supplied-tree control, with candidate groups, membership constraints,
raw gene trees and reconciliation rules held fixed. The control exactly
reproduced the inferred baseline partition and74.106074% F1. Three variants
had rooted clade distance2 and three distance4. Their F1 ranged from73.744102%
to74.270946%, or -0.361972 to+0.164872 percentage points relative to control.
Every one of the18 Bonferroni-adjusted F1/precision/recall intervals included
zero. Family F1 tied in62-69 of70 RefOGs per perturbation, although complete
predicted partitions changed. All runs preserved251,378 genes and passed
native-provenance checks and official-score crosschecks.

These observations describe this limited, development-exposed topology panel;
they do not establish equivalence, robustness to arbitrary species-tree error,
or improved tree estimation. Intervals remain compatible with larger effects
than the observed differences, and the variants are not posterior draws.
No variant was selected to change the method. Recorded reconciliation costs
reuse raw gene trees and are not end-to-end efficiency measurements.
[Full scores, uncertainty and coverage](OB_SPECIES_TREE_ROBUSTNESS_RESULTS_20260916.md).

![Supplied-tree perturbation scores and all paired effects](figures_species_tree_robustness_20260916/species_tree_robustness.png)

Figure: control and six fixed topology perturbations, ordered by the frozen
variant indices. Thick intervals are nominal95% paired RefOG intervals;
thin intervals retain Bonferroni adjustment over all18 endpoints using20,000
paired resamples. Gene-tree checkpoints are shared; rooting/reconciliation
is recomputed against each supplied species tree.

### Exploratory Feature Strata Show Precision-Recall Tradeoffs

Fourteen feature strata were frozen before joining the new features to method
errors. They describe reference-family size, copy number, canonical alignment
identity, relative sequence length and global composition, not independently
verified duplication history, fragments or domain architecture. All70 RefOGs
remain represented in each dimension. Twelve families have missing identity
because at least one aligned pair lacks comparable canonical positions.

Both OrthoHMM configurations were compared with full OrthoFinder3.1.5 for
F1, precision and recall, retaining84 endpoints in multiplicity adjustment.
Eleven bins had at least five families, yielding66 interval-bearing endpoints;
one single-family composition bin was descriptive only and two empty bins
were explicitly nonestimable. None of the22 adjusted F1 intervals excluded
zero. Nine adjusted precision/recall intervals excluded zero: two precision
advantages and seven recall deficits for OrthoHMM across overlapping strata.
These are not independent replications or tests of differences between strata.

In the40-family short-relative stratum, high-sensitivity recall was18.901
percentage points below OrthoFinder (adjusted interval[-36.475,-5.208]);
phylogenetic-mode recall was14.833 points lower ([-34.240,-0.370]). In the
six large families, observed F1 differences were-16.469 and-11.616 points,
but both adjusted F1 intervals included zero. These observations do not prove
a fragmentation mechanism or justify outcome-selected parameter changes.

All three native-format parsers preserved251,378 input genes, and fresh
full-reference sufficient statistics and official-score crosschecks reproduced
the frozen overall estimates. Stratified estimates restrict those original
per-family counts without changing the reference or low-certainty conventions.
Percentile intervals are approximate; the84-endpoint adjusted tails contain
only about six of20,000 resamples each. [All effects and limitations](OB_STRATIFIED_ERROR_RESULTS_20260916.md).

![All stratified OrthoBench F1, precision and recall effects](figures_ob_stratified_errors_20260916/stratified_errors.png)

Figure: both OrthoHMM configurations minus full OrthoFinder3.1.5 across all14
prespecified strata, with identical horizontal scales. Filled points show66
interval-bearing endpoints: thick lines are nominal95% paired intervals and
thin lines retain the84-endpoint Bonferroni adjustment. Six open points are
descriptive estimates from a single-family stratum;12 endpoints in empty strata
remain explicitly nonestimable. Family counts are shown for every row. Strata
overlap, and differences between their effect estimates are not interaction
tests. Scores and points are checked against full-reference sufficient statistics;
the figure does not add new inference runs or independent biological evidence.

### Retained Checkpoints Localize Membership Changes

All70 OrthoBench families were traced across six retained checkpoints, with
40,733 pair rows checked against native partitions. Both assignments of the
one shared reference gene were preserved. The8,440 logged candidate merges
reconstructed the candidate partition;168 events touched reference genes.
These checks establish consistency of the retained evidence, not biological
correctness. Profile expansion starts from unrefined multipass clusters;
profile-off refinement is a separate branch, not its input.

Comparing matched refined endpoints, the profile-on branch gained836 and
lost160 within-family pairs. Candidate expansion added4,272; final tree,
reconciliation and constraint processing removed1,575. These raw co-membership
counts include within-species pairs and low-certainty assignments and are not
official benchmark recall. All six hash-selected illustrations were retained,
including neutral examples and candidate gains subsequently partly removed.
[Trace interpretation and full-family evidence](OB_FAMILY_TRACE_INTERPRETATION_20260916.md).
Rejected-search candidates, added profile-edge identities and per-family
tree/constraint mechanisms still require inspection; no biological mechanism
or new default is inferred from these membership changes alone.

Reconstruction of the250 reference-incident candidate families distinguished
root-lineage extraction from subsequent satellite membership constraints.
For114 reconciled and136 bypassed candidates, saved node calls and all193
local constraints exactly reproduced native final groups. Of1,575 lost
within-reference pairs,134 were separated by root-lineage rules and1,441
by subsequent constraints, affecting four and six families respectively
(one family in both sets). Among the frozen illustrations, losses inRefOG014
andRefOG021 occurred at the constraint step. These are execution-order
descriptions, not official recall, independent tree validation or a claim that
constraints should be removed; the unconstrained control showed a precision
cost without established F1 improvement.
[Full reconciliation trace and limitations](OB_RECONCILIATION_TRACE_20260916.md).

### Limited Parameter Sensitivity

The prespecified OrthoBench neighborhood changed CPM resolution, candidate
minimum normalized score, or candidate minimum margin by20% one at a time,
retaining all other settings and the full inferred-tree pipeline. CPM variants
recomputed their own HMM-refined seeds and candidates; threshold variants reused
the verified baseline seeds. Exact-input raw-tree reuse was allowed, but supplied
species trees were not substituted. All seven configurations retained251,378
input genes, passed native admission, and matched the official OrthoBench scorer.

Control F1 was74.106074%. F1 values were74.586147% and71.463468% for CPM0.08
and0.12,74.008570% and74.110826% for minimum normalized scores0.024 and0.036,
and74.973144% and73.274559% for minimum margins1.2 and1.8. Each of the six
F1 contrasts had a Bonferroni-adjusted interval including zero, using20,000
paired RefOG draws and all18 planned F1/precision/recall endpoints. The only
adjusted interval excluding zero was the recall decrease at CPM0.12:
-2.010 percentage points, adjusted interval[-6.291,-0.148]. Minimum normalized
score changes affected only one or three families' F1 values and left aggregate
recall unchanged; this is a local observation, not general invariance.

The largest observed F1 increase,0.867 points at minimum margin1.2, had an
adjusted interval[-1.362,4.316]. It does not justify default promotion or
an OrthoFinder superiority claim. This post-development panel is not independent
confirmation, and intervals including zero do not establish equivalence.
Shared-node phylogeny timings are incremental and omit upstream replay/candidate
costs. The corresponding QfO panel awaits a reproducible native baseline.
[Complete scores and paired intervals](OB_PARAMETER_NEIGHBORHOOD_RESULTS_20260916.md).
[Parameter sensitivity figure](figures_ob_parameter_neighborhood_20260916/parameter_neighborhood.pdf).

## Limitations And Unfinished Analyses

### Remaining Requirements

No universal superiority, arbitrary-dataset generalization, or controlled
speedup is established. Frozen YGOB evaluation is complete and supports
novel-taxon transfer, not family-disjoint confirmation. Homolog-family
overlap remains substantial even after taxon exclusions. The proposed HMM
contribution is not established by the completed exploratory OrthoBench
sequence-search control; QfO extension remains outstanding. The interaction
between broader candidates and reconciliation has descriptive OrthoBench
evidence but awaits QfO evaluation and additional controls. Corrected multi-seed simulations are complete
but do not establish an OrthoHMM advantage or profile-expansion benefit.
QfO and independently annotated error strata, mechanistic tracing, broader tree-error
robustness and the QfO parameter panel, matched resource scaling,
and a prespecified biological application remain required.

QfO historical replay equivalence remains unresolved. A bounded capture
verified identical initial RBNH graph arrays and gene order but a different
first clustering partition than an earlier diagnostic. Thus the observed
discrepancy precedes profile refinement; its specific cause is not yet
established. No replay result is silently substituted for the historical
baseline. [Diagnostic evidence](QFO_REPLAY_DRIFT_DIAGNOSIS_20260916.md).
Three subsequent instrumented single-CPU workers with matching recorded
native-library/runtime identities reproduced the earlier diagnostic partition,
but not the replay capture. A subsequent controlled affinity panel exposed
within-configuration variability: the two one-CPU repeats differed (349,898
versus349,950 groups), although their recorded sources, libraries, inputs,
settings and affinity matched. Both32-CPU repeats produced349,898 groups.
Independent file and complete-partition admission confirmed the discrepancy.
Thus CPU availability alone does not explain the variation, and the earlier
matching repeats do not establish general repeatability. The cause remains
unresolved; no preferred partition is substituted for the historical baseline.

A subsequent native-boundary diagnostic stopped during its second worker
before optimization: ordered native endpoint hashes differed from the saved
graph despite matching vertex/edge counts and ordered weights. The first worker
had matched the saved graph before and after optimization. Fresh saved-array
hashing confirmed the expected endpoint fingerprint. This localizes a discrepancy
to pre-optimizer graph handling in that instrumented run, but does not yet
distinguish conversion, construction or observation effects, or establish the
cause of all historical partition differences. The failed partial panel is
preserved and not presented as a completed reproducibility experiment.

A second boundary diagnostic recorded six incorrect native edges despite an
intact original constructor array. A completed construction-only follow-up
alternated three original-int32 workers with three explicit-int64-copy workers.
Five graphs matched the saved endpoint stream; one explicit-int64 worker had
the same six endpoint mismatches. Both its original and converted arrays were
intact. Native tuple and source/target access agreed on the observed endpoints,
and edge lookup did not find the expected pairs. Independent admission checked
257 file records and reconstructed each worker's full reported native endpoint
hash from the saved arrays plus its complete bounded mismatch witnesses.
Explicit conversion alone is therefore not a sufficient fix in this experiment.
The audit validates preserved observations, not an independent inspection of
the historical live object; it does not establish a specific library defect,
hardware cause, or the cause of all historical variability. No optimizer was
invoked and no partition or accuracy score was produced in this follow-up.
[Construction evidence and limitations](QFO_CONSTRUCTION_DIAGNOSTIC_20260916.md).

A direct-construction follow-up completed six workers alternating minimal and
frozen-worker imports. Four had the same six endpoint mismatches before weight
assignment, two matched, and each worker's recorded differences were unchanged
after assigning weights. Both modes contained mismatching workers with intact
constructor arrays. Independent admission checked278 file records and reconstructed
all post-weight endpoint hashes. Pre-weight hashes are implied by bounded complete
witnesses, not separately recorded full hashes. Neither weight assignment nor
OrthoHMM/Leiden imports is necessary for this observed failure; its underlying
cause remains unresolved. This is preserved-observation evidence, not independent
live-object inspection or proof that every historical result is affected.
[Direct-stage results](QFO_DIRECT_GRAPH_RESULTS_20260916.md).

Historical timing and memory records differ in scope and accounting.
Cached replays are incremental computations, not end-to-end timings;
supplied trees have upstream costs. Unscheduled competing workloads were
observed, so exclusive Slurm allocation alone cannot establish controlled
CPU conditions. GNU-time maximum RSS and simultaneous process-tree RSS
must not be pooled as if they were the same measurement.

## Data And Code Availability

Code and analysis milestones are tracked in the authorized repository.
Machine-readable reports retain input, source, and output hashes where
audited; large raw datasets and working outputs are not included in normal
source commits. The figure bundle is generated from the comparison JSON
and includes a provenance manifest. Absolute local paths are provenance,
not portable download locations. A portable workflow bundle, complete
dependency lock, raw-output inventory, redistribution-permission review,
versioned release, and archival deposition remain unfinished. No archive
accession or publication DOI has been assigned by this work.

See the [claim-to-evidence checklist](PUBLICATION_CLAIMS_20260916.md) before
strengthening conclusions or describing this package as publication-ready.
