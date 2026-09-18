# Recovered QfO Factorial Protocol

This development-exposed experiment extends the completed OrthoBench factorial
to QfO. It is frozen before new candidate expansion or reconciliation results
are inspected. Historical QfO scores and the four recovered-stage scores are
already known; this is not prospective independent validation.

## Design

Use all 976,504 genes in the 78 validated QfO proteomes and the admitted
recovered v2 normalized-HMM-hit checkpoint. The frozen core remains 7f3a9e4,
with its corrected installed native runtime and source-equivalent isolated
launcher. Do not substitute the later development banding correction.

Eight cells cross three binary factors:

- P: profile branch off/on, starting respectively at recovered
  `multipass_refined` and `strict_profiles_refined` complete partitions.
- C: no additional candidate expansion / frozen `satellite_v2` expansion.
- R: complete candidate partition / inferred-tree reconciliation with
  `species_overlap`, `positive_paralogy`, minimum-variance species-tree rooting,
  MAFFT and FastTree. C-on/R-on retains native seed membership constraints.

P-off still uses the initial HMM search. The profile contrast includes
downstream sequence refinement; it is not an isolated change in profile scores.
C changes both candidate membership and, when R is on, applicable seed-based
constraints. R-on includes species-tree inference from its own candidate arm;
do not label this a fixed-tree-only reconciliation contrast. R-off predictions
are cross-species clique pairs from complete candidate groups; R-on conversion
must follow the already audited native RootHOG semantics. Neither may be
silently replaced with pre-clustering edges or hierarchical-group exports.

Preparation uses the existing OrthoBench candidate and cell-planning helpers,
but QfO has a separate admission gate. The recovered replay is NOT asserted
equivalent to historical comparison rows. Check all retained replay provenance,
the exact numeric checkpoint, FASTA ownership, complete partition coverage,
and membership-constraint validity before launching reconciliation. Preserve
failed preparations and never overwrite completed arms. Record source/runtime
hashes before and after processing. Freeze and verify the reconciliation
launcher and its dependencies before execution; generated cell commands alone
are not admission evidence.

## Endpoints and Contrasts

Retain all six native QfO endpoints, individual precision/recall where defined,
scored relation coverage and the explicitly secondary project-defined mean.
Reuse existing recovered C-off/R-off scores only after exact partition,
conversion, reference and scoring provenance equality is demonstrated. Do not
transfer scores from historical high-sensitivity or satellite rows.

For SwissTrees, primary paired inference uses the same 18 native families,
raw/2+1 native confusion-count rule, family-macro precision and recall, and
harmonic F1. Recompute those aggregates in each of 100,000 shared family
bootstrap draws, NumPy PCG64 seed 20260922. Retain all twelve one-factor
simple effects (on minus off at every setting of the other two factors),
plus C-by-R differences of differences separately for each P level. For each
of these fourteen contrasts report F1, precision and recall: 42 endpoints,
nominal 95% percentile intervals and Bonferroni-adjusted percentile intervals
at .05/(2*42) and 1-.05/(2*42), linear quantiles. Retain family-level paired
differences and wins/ties/losses. Conditional development-family intervals are
not selection-adjusted or proof of generalization. No extra contrasts will be
added after outcomes without explicitly labeling them exploratory.

Other QfO challenges retain native descriptive outputs. Appropriate independent
units for their paired uncertainty remain unresolved; no IID pair bootstrap
or pooled use of unlike native stderr fields is permitted. The six-metric
mean has no implied valid joint confidence interval. No tuning or default
promotion is authorized by this experiment. Negative/neutral outcomes and
failed cells remain visible; failure is not assigned an accuracy of zero.

## Resources and Execution

Prepare on `bizon`, 32 CPUs/192 GiB, with hash seed zero and single-threaded
BLAS/OpenMP. This does not use the exclusive DGX timing node. Record incremental
preparation and reconciliation costs, scheduler outcomes and resource evidence,
but do not pool shared-workstation replay costs with controlled end-to-end
DGX timings. All four reconciliation cells require separate native-output
and input/runtime admission before conversion/scoring. Preparation alone does
not complete the ablation or establish any accuracy claim.
