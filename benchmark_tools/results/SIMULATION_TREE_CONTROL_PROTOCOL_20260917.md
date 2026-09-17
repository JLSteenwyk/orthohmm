# Simulation Generating-Tree Controls

This follow-up addresses the outstanding simulation-truth species-tree-error
requirement. Previous variable-length outcomes have been inspected; this is
development-exposed mechanism/robustness analysis, not independent confirmation
or a search for defaults that favor OrthoHMM.

## Fixed Panel

Retain all70 cells of the frozen variable-length panel: ten seeds20261101-20261110
and baseline, divergent, turnover, divergent_turnover, missing20, uneven_taxa,
taxon_count_control. Do not select successful original method runs or replace
failed seeds. Verify generation/history/input manifests and original sequence
bytes before using any supplied tree. Source manifest SHA256:
806aa1e5f6976c323ff2f6641dd88e25264e7417749e7d77565294666aee776b.

For each dataset, take its parent's simulated `T/ExtantTree.nwk` and prune it to
the exact retained FASTA species basenames. Require unique taxa, bifurcating
rooted topology, nonnegative finite lengths, and unchanged retained patristic
distances after pruning. If the old root becomes unary, use the induced rooted
subtree. No protein sequence, gene-family truth or model parameter is modified.

Prepare three trees:

1. The induced generating tree, a supplied oracle diagnostic.
2. The first SHA256-ranked one-NNI neighbor, rooted clade distance2.
3. The first SHA256-ranked distance4 topology among neighbors of the first three
   SHA256-ranked one-NNI trees, using the existing OrthoBench perturbation helper.

Record helper/source/tree hashes, topology digests, parent topology and every
taxon. Remove internal labels/support from serialized supplied trees. Lengths
travel with moved subtrees: perturbed trees need not remain ultrametric, preserve
patristic distances, or represent an empirical error distribution. The two-NNI
tree need not descend from the single one-NNI tree retained as an analysis arm.

## Execution and Interpretation

Evaluate frozen satellite_v2 and full OrthoFinder3.1.5 with each supplied tree:
420 planned method/tree/cell combinations, in addition to their existing
inferred-tree baselines. Keep high-sensitivity results as a descriptive
phylogeny-free reference, not another tree-perturbation arm. Inapplicability and
failures must remain in the complete inventory rather than disappearing from
the design. Known pre-phylogeny comparator failures may persist with supplied
trees and do not demonstrate a tree-inference problem.

Before execution, pin exact supported CLI commands and output conversions.
Validate supplied-tree naming/root semantics and unchanged input/candidate
coverage. Reuse expensive upstream checkpoints only after verifying their
identity and that reconciliation is actually recomputed. An unchanged inferred
tree supplied through the new path must reproduce its original downstream
partition before attributing effects to topology. If this mode gate fails,
investigate it without tuning or choosing favorable trees. Failed original
inference must not be excluded merely because no inferred-tree control exists.

Generating trees reveal information unavailable to ordinary inference. They
are not an achievable end-to-end baseline or necessarily an accuracy upper
bound in the presence of gene-tree/model errors. Separate changes in tree
availability, upstream grouping, rooting and reconciliation where evidence
allows. Do not assume any observed gain identifies the cause of prior failures.

## Analysis

Preserve the existing simulation truth definition and native-pair/output
conversions. For each condition and method, compare generating-tree versus
inferred-tree F1/P/R, and each perturbation versus generating-tree F1/P/R.
This specifies7 conditions x2 methods x3 contrasts x3 endpoints =126 exploratory
endpoints. Keep this multiplicity count even when cells cannot be evaluated.

Use paired resampling of the ten simulation seeds, not individual dependent
gene pairs, with20,000 replicates and fixed seed20260918. The statistic is the
arithmetic mean of per-dataset metric differences over paired successful seeds;
recompute that mean in each replicate, without pooling gene pairs across datasets.
Use two-sided Bonferroni intervals with family alpha0.05 over126 endpoints.
Report paired complete-case
sample counts and all failures separately. Report effects, intervals and seed-
level results without selecting a best topology or changing defaults. These
small simplified panels do not validate arbitrary species-tree error robustness,
real biological domains/fragments, or independent generalization.

Tree preparation alone is not completed inference or a scored experiment.

## Prepared Inputs

All70 datasets passed generation/history/input checks and tree preparation, with
no inapplicable cell. There are210 trees and420 unique planned method/tree runs.
Independent rereading verified each tree hash, exact retained taxa and rooted
clade distance against the generating tree projected onto those taxa, separately
from the pruning implementation. The
[prepared manifest](simulation_tree_controls_prepared_20260917.json) has SHA256
0caffab16f73019fabafe1fcfaac8c2de90960c5f1317803bf83308a6d2bf914.
No supplied-tree inference or scoring has been launched. CLI/output semantics,
mode-equivalence checks and native execution still precede scientific admission.
