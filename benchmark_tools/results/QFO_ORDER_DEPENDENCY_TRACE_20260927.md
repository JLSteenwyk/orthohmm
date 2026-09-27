# QfO Ordering Dependency Trace

The [machine-readable trace](qfo_order_dependency_trace_20260927.json)
uses the installed, pinned phylogeny implementation and checks all 78 input
FASTA identities and candidate-stage receipts before and after inspection.
It runs selection only, not alignment, tree inference, reconciliation or scoring.

- All 351,737 shared candidate families retain their original family IDs.
- `Family0014602` changes from a singleton to a 48-gene family.
- `Family0015042` changes from a 48-gene family to a singleton.
- The singleton changes from `sp|Q8TAF7|ZN461_HUMAN` to
  `tr|A0A2I3RT74|A0A2I3RT74_PANTR`; the multi-gene family swaps those members.
- Each arm has one reconciliation-requiring changed family and one bypassed
  singleton. The multi-gene family covers 16 species.
- Five constraints differ on each side as a multiset, so the earlier five
  unequal positions are not merely a permutation of identical constraints.
- Both arms select the exact same ordered 26 marker families and gene lists.

This bounds the changed direct inputs but does not establish species-tree,
gene-tree, root-group, native-pair, or benchmark-score equivalence. Reusing
historical outputs without validating dependencies would still be unjustified.
No reference labels were used. Thirty-four focused trace/readback/launcher
tests pass. No scientific default or historical score changed.

## Frozen Downstream Comparison

Before any downstream native execution, fix the following scientific design:

1. Run the installed recovery phylogeny stage on the retained-order candidate
   partition with its constraints, using all 78 admitted proteomes. Infer
   species and gene trees afresh, without historical tree checkpoints.
2. After structural and artifact validation of that arm, run the same stage
   on canonical candidates and their constraints. Reuse only copied raw-tree
   checkpoints and associated alignments/FASTA artifacts whose sequence,
   membership, tool bytes, configuration and checkpoint hashes validate.
   Recompute rooting, reconciliation and membership filtering for every family.
   Infer the species tree afresh in the second arm; do not copy its cache.
3. Use `species_overlap` root rule, `positive_paralogy` pair rule,
   `min_variance` inferred species-tree rooting, the same validated recovery
   environment, MAFFT 7.525 and the retained FastTree 2.2 AVX2 executable.
   Pin complete executable/helper identities, input order, sources, arguments
   and environment in an execution plan before submission.
4. Request 32 local CPUs and 128 GiB, at most six hours per arm, one attempt
   each with no automatic retry. No search rerun or parameter optimization.
5. Validate full input coverage, structural output consistency and provenance.
   Compare both arms and historical outputs in all three contrasts: complete
   root partitions, native ortholog-pair sets, pair confidence, inferred species
   tree, reconciliation counts and membership audit. Report family-level
   differences and checkpoint reuse, not only aggregate counts.
6. Preserve every difference. Root/pair equality does not prove equality of
   other products. Do not transfer historical QfO metrics when their actual
   evaluated prediction sets differ; freeze and execute the corresponding
   existing scoring workflow separately if rescoring is needed.

This is a downstream cached-search reproducibility comparison, not a fresh
end-to-end search validation or independent generalization experiment. Resource
measurements on the shared host are descriptive only. Execution preparation,
native runs and readback remain unfinished at this milestone. The full
publication requirements and DGX deferral remain unchanged.
