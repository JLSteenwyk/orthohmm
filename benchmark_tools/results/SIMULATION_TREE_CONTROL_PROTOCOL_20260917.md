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
# Native Input Compatibility Audit (2026-09-17)

The independently admitted unchanged-tree pilot now supports expansion to the
full mode-control panel. Frozen inventory
simulation_mode_panel_prepared_20260917.json SHA256
03f0d115d0f2a1f5354ca6a8a369e419dbd9e10945e450dfd2534d61212bb747
contains all70datasets/140method slots:130new method runs,2reused pilot runs,
and8unavailable inferred-baseline controls (3OrthoHMM,5OrthoFinder). Availability
uses original completion status only, not scores. The other method still runs
when one baseline is unavailable. These8failures remain in the420oracle/topology
run inventory; missing mode controls do not authorize a tree-only causal claim.
The array executes at most2four-CPU16GiB tasks at once on the shared host. These
are scientific controls, not matched scaling measurements. Every task keeps
unavailable-method records and records new failures; the whole panel needs
independent native validation before the main tree-effect contrasts are scored.

Before execution, direct checks with the installed native parsers found that
OrthoFinder 3.1.5 reads the original leading `[&R]` serialization as a `NoName`
tree, rejecting the expected taxa. Frozen OrthoHMM accepts all 210 originals.
Do not pass these original files directly to OrthoFinder, and do not overwrite
them or change the frozen preparation manifest.

`prepare_portable_simulation_trees.py` creates derivative plain-Newick files
using Bio.Phylo, without the leading annotation. The binary root is retained;
every descendant clade and its branch length must match exactly after rereading.
Both OrthoFinder's `CheckUserSpeciesTree` and frozen OrthoHMM's
`parse_species_tree` accepted all 210 derivative files in direct local checks.
This is parser acceptance, not reconciliation or prediction equivalence.

Derivative manifest: `simulation_portable_trees_prepared_20260917.json`, SHA256
`b9ed4fb8dc27da28dd56c674d1ece2edbb3a04697dec14ea3d4538bf6d9dbc0b`.
Files remain at `benchmarks/results/simulation_portable_trees_v1`.
The machine-readable preparation status intentionally does not claim native
inference validation. Recheck parsers and all source/tree hashes in execution
preflight; the direct parser checks above are recorded here, not an independent
admission artifact.

The tested `simulation_supplied_commands.fresh_supplied_method` constructs
isolated fresh runs from the frozen method command dictionaries. OrthoHMM changes
only the output locations, tree mode to `supplied`, and `--species-tree` path.
OrthoFinder changes the copied FASTA location and adds `-s`; all other baseline
arguments remain unchanged. Existing outputs, overlapping input/output paths,
preexisting supplied trees and restart flags are rejected. The full 420-row
inventory was dry-constructed with unique destinations; none was executed.
Use the portable derivative tree paths when materializing actual commands.

Local source audit used the installed OrthoFinder 3.1.5 package:
`run/species_info.py:282` checks unique exact taxa and a binary root;
`comparative_genomics/orthologues.py:697` (`RootSpeciesTree`) bypasses STRIDE
and its multiple-root handling for supplied trees. Frozen OrthoHMM
`orthohmm/phylogeny_pipeline.py:1050` parses and writes a supplied tree instead
of calling its tree-inference branch. Consequently, unchanged supplied-tree
controls remain mandatory: mode changes can affect more than serialization.
The `-ft` restart route has not been validated for this panel; do not assume
postprocessed result directories recover the original upstream inference state.
Fresh-run candidate and gene-tree identity must also be checked before treating
differences as tree-only effects. No equivalence gate has been waived.

## Main Experiment Runner

`run_simulation_tree_experiment.py` and the prepared210task batch implement the
420method runs (two tools per supplied-tree cell). They are not submitted until
mode-panel admission21367 finishes and its exact result hash is reviewed and
supplied. The runner requires all132available controls to be independently
equivalent and all8unavailable controls to match their original failure records.
Any new failed/non-equivalent mode control stops authorization for review; it
does not drop that dataset, select a different seed, or relax validation.

Every cell rechecks the frozen original/portable tree identity, exact taxa,
topology and branch lengths, generation and sequence inputs, both native tree
parsers, tool environments and native outputs. Both methods execute even when
one had no original inferred-tree baseline. Fresh outputs retain method failures
and whether the requested rooted topology reached native results. No truth
scoring or parameter selection occurs inside this runner. Retained candidate,
graph, gene-tree and alignment records support the later cross-arm audit.

Dry construction checked all210portable tree cells and420unique fresh commands.
No main experiment inference has started. The batch uses at most2concurrent
4CPU16GiB tasks, writing `simulation_tree_experiments_v1/cell_<index>`. Input
preparation and native-run validation are not independent scientific admission;
the latter and paired126-endpoint analysis remain separate required steps.

## Complete-Panel Admission (17 September Update)

The earlier preparation-only status above has been superseded: all available
unchanged-tree controls passed independent pair and downstream partition checks;
see SIMULATION_MODE_PANEL_VERIFIED_20260917.md. Main cell0 ran as21405_0 and
cells1-209 are executing as21406. Neither array is a matched timing experiment.

`admit_simulation_tree_panel.py` requires all210 scheduler tasks uniquely terminal
across those exact two arrays before it creates an admission report. Missing,
duplicated, nonterminal or contradictory completion records reject. Unsuccessful
terminal tasks remain in the420-method inventory rather than being dropped.

For completed cells the auditor reconstructs non-tree commands independently of
the execution command builder, checks source/helper revisions, frozen manifests,
mode-control evidence, tree identity and lengths, and reruns input, environment,
native-parser and native-output validation. Native prediction adapters verify
their inputs against the recorded output inventory. Root topology and retained
upstream artifact inventories are independently recomputed and compared to the
runner's observations. Metadata disagreement rejects admission rather than
silently reclassifying a cell. A genuinely unretained requested topology remains
an explicit outcome, not a tree-accuracy estimate.

This admission does not score evolutionary truth. Cross-arm upstream artifact
comparisons and the prespecified126-endpoint analysis remain separate steps;
valid native execution alone does not establish a tree-only causal effect.
