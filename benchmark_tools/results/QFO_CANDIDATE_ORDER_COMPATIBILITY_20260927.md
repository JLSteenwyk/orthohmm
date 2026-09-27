# QfO Candidate-Input Ordering Compatibility

The [read-only audit](qfo_candidate_order_compatibility_20260927.json) binds
the retained corrected-QfO candidate admission, preparation manifest and
numeric checkpoint, checking all checkpoint file hashes before and after
inspection. It uses memory-mapped arrays and bounded chunks; it does not
change arrays, predictions, parameters or scores, and launches no inference.

The 984,137 gene identifiers are unique and lexically sorted. All 90,687,327
hit rows have valid indices and finite positive scores, including 983,835
self-hits. Their query/target order contains **34,986,309 adjacent descents**
and zero adjacent duplicate pairs. Because the order is not sorted, zero
adjacent duplicates does not prove global uniqueness.

Thus the new canonical directed-pair policy is not an identity transformation
of this retained input. This does **not** establish changed candidate groups,
changed reconciliation or changed accuracy. Nor does it establish failure of
the retained historical results. Those remain unchanged and must not be
relabeled as validation of the new entrypoint.

The original admitted corrected-QfO high-sensitivity job 21707 took 19:41:39
on 32 allocated CPUs. Rather than repeat that expensive search merely to test
ordering, the next experiment should compare historical and canonical order
using the same admitted hits and seed groups, frozen expansion implementation
and fixed runtime. Record whole candidate partitions and semantic membership
constraints before considering fresh phylogeny or rescoring. Preserve any
differences and do not tune to recover historical scores. This candidate-only
experiment will not replace eventual installed end-to-end QfO validation.

Fourteen focused tests pass, including chunk-boundary descents/duplicates,
within-query target ordering, nonadjacent-duplicate limitations, invalid scores,
empty arrays and comparison with an independent Python pair-order calculation.
The publication goal remains incomplete; no new accuracy claim is made.
