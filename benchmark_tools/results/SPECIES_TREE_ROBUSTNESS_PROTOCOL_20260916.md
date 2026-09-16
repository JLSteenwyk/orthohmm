# Species-Tree Robustness Diagnostic

This is development-exposed robustness analysis, specified after earlier
OrthoBench/YGOB results. It does not tune the frozen method, establish
independent confirmation, or use YGOB outcomes to select tree changes.

## Fixed OrthoBench Panel

Use the validated p1_c1_r1 rooted species tree and the same candidate groups,
membership constraints, sequence inputs, gene-tree tools and reconciliation
rules. Supply seven explicit rooted trees: an unchanged-topology control,
three rooted nearest-neighbor interchanges (NNIs), and three two-NNI variants.
The latter have rooted clade symmetric-difference distance four from baseline;
one-NNI variants have distance two. Retain every taxon exactly once.

Enumerate rooted binary NNI neighbors; choose the first three unique
topologies sorted by SHA256 of canonical nontrivial rooted clades. Enumerate
second neighbors from those three, restrict distance to four, and select
three by the same rule. No reference-family outcomes enter tree selection.
Record exact input/output hashes and the parent topology for every variant.

Lengths travel with their subtrees. Remove internal confidence labels, which
do not measure support for the perturbed trees. These edits are controlled
topological stress tests, not draws from a posterior or an empirical error
distribution. They do not preserve all evolutionary distances.

## Execution And Admission

Use the frozen scientific core/native runtime and supplied-tree mode. Preserve
species_overlap root rule, positive_paralogy pair rule and high-confidence
membership constraints. Give every run a new output directory. Verified raw
alignment/gene-tree checkpoints may be reused, but reconciliation must be
recomputed against the supplied tree; validate this behavior before launching.

First require the unchanged supplied-tree control to reproduce the inferred
p1_c1_r1 root-HOG partition. If it does not, investigate mode conversion or
cache behavior before attributing changes to topology. Do not overwrite the
original inferred-tree output or call replay costs end-to-end timings.

After native/provenance/taxon/partition checks, use official OrthoBench scoring
with paired70-RefOG bootstrap20,000 replicates, seed20260918. Compare each of
six perturbations to supplied control for F1/P/R:18 exploratory endpoints,
retained together in Bonferroni adjustment even if a perturbation fails.
Report failures, all effects and coverage; no selection of a best tree.

This panel does not satisfy the separate parameter-neighborhood, QfO
robustness, simulation-truth tree-error or controlled scaling requirements.
Tree generation alone is not a completed reconciliation/accuracy experiment.
