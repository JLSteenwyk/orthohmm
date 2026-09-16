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
multiplicity handling. No accuracy outcome is available in this draft.
[Frozen protocol](YGOB_VALIDATION_PROTOCOL_20260916.md).

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
These experiments are specified but not completed.
[Ablation protocol](PUBLICATION_ABLATION_PROTOCOL_20260916.md).

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

## Limitations And Unfinished Analyses

No universal superiority, arbitrary-dataset generalization, or controlled
speedup is established. Independent validation is pending. Homolog-family
overlap remains substantial even after taxon exclusions. The proposed HMM
contribution requires a matched sequence-search control, and the interaction
between broader candidates and reconciliation remains unmeasured under the
prospective factorial design. Error strata, mechanistic tracing, validated
multi-seed simulation, tree/parameter robustness, matched resource scaling,
and a prespecified biological application remain required.

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
