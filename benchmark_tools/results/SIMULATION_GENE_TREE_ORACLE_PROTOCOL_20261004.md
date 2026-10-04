# Fixed-Candidate Generating Gene-Tree Diagnostic

## Scope Frozen Before New Outcomes

Address goal4.4: the retained species-tree experiment cannot distinguish
gene-tree estimation/rooting from reconciliation errors. Reuse every one of
the70 generating-species-tree OrthoHMM satellite_v2 records (seven frozen
conditions, ten seeds20261101-20261110) in the independently admitted native
panel. Do not select datasets or candidates by accuracy. Prior simulation
outcomes are development-exposed; this is an exploratory mechanism diagnostic,
not independent confirmation, a new default, or an end-to-end benchmark.

Pin admission17a1be71823b9fb7fb13082d05001d7667bb89c3366da2e3138a124e90907820,
tree preparation0caffab16f73019fabafe1fcfaac8c2de90960c5f1317803bf83308a6d2bf914,
and the frozen phylogeny source216d97608e6dede8960f3b06fe7036bd23da8fae5677df88583ba9369ec3c1bf.
Check selected native reports, execution inventories, candidate partitions,
gene/species trees, checkpoints, native pairs, constraints and truth inputs
against retained identities. Generating gene trees must match the original
parent truth inventory, including for derived missing/taxon datasets.

## Fixed Intervention

Hold the supplied species tree, candidate membership, reconciliation rules
(species_overlap/positive_paralogy), and satellite-constraint policy fixed.
Recover owners from pinned FASTA files and ancestral-family membership from
the retained truth, not heuristic parsing of inferred names.

Every candidate remains in the inventory. Only a candidate that actually
underwent gene-tree inference and belongs to exactly one generating ancestral
family is eligible. Mixed-family candidates have no single generating tree:
retain native predictions, mark ineligible, and never split them using truth.
Unambiguous bypass candidates retain native predictions without intervention.
Any unavailable native cell or failed control is an explicit outcome, not an
excuse to choose another seed or silently drop a candidate.

Recompute the frozen baseline from each retained rooted gene tree and require
exact native pair-set reproduction after the original constraint policy. For
eligible candidates induce the rooted generating tree on the exact candidate
genes with DendroPy, suppressing unary paths and preserving the induced root.
Rename only extant leaves to their validated global FASTA identifiers; remove
ancestral internal labels (not evidence of inferred branch support).

Compare three fixed arms:

1. Retained inferred rooted gene trees (baseline control).
2. Induced generating gene tree with its generating root.
3. The same induced generating tree rerooted by frozen minimum-duplication/loss.

Reconciliation and constraints are recomputed; bypass and ineligible predictions
are unchanged. Record native-pair additions/removals and TP/FP/FN/F1/P/R per
dataset and candidate. Use the existing cross-species simulation truth and
pair scorer; do not replace pair predictions with cluster-derived pairs.
Record rooted clade and unrooted split distances between inferred and induced
generating trees using DendroPy's shared taxon namespace. Zero unrooted distance
with nonzero rooted distance identifies a root disagreement, not necessarily
its causal effect on the endpoint.

## Reporting And Limits

Report all70 cells, candidate eligibility and exact finite-panel per-condition
means; do not claim population intervals or significance from descriptive
changes. No best arm is selected. Keep gains, losses and unchanged outcomes.
True trees and roots are unavailable in real inference. The intervention also
changes branch lengths/support information, so it is not a topology-only causal
experiment or an accuracy upper bound. Pair endpoint improvements do not prove
correct ancestral-copy/root-HOG membership. Residual oracle errors can reflect
loss-obscured duplications, reconciliation assumptions, upstream candidates and
constraints; do not infer their biological cause from F1 alone.

Use DendroPy5.0.8/Biopython1.87 in the retained local reporting environment,
one numerical thread. This is a small read-only diagnostic on existing inputs,
not new HMM search, native inference or matched timing. Do not disturb other
analyses or modify existing results. Preserve input hashes at completion and
write a new result only. Source/protocol commits precede actual execution.

API semantics reviewed against the local5.0.8 implementation and official
[tree manipulation](https://jeetsukumaran.github.io/DendroPy/primer/treemanips.html)
and [tree comparison](https://jeetsukumaran.github.io/DendroPy/library/treecompare.html)
documentation. These support tree operations, not the scientific accuracy of
OrthoHMM or the proposed intervention.
