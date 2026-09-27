# Paired QfO Candidate Ordering Result

Job 22328 completed successfully in 00:03:08 on two local CPUs. Both arms
used the same admitted 90,687,327 hits, 984,137 genes and seed partition,
installed implementation, environment and candidate settings. Only hit order
changed. The [protocol](QFO_ORDER_REPLAY_PROTOCOL_20260927.md) was committed
and pushed before submission. The [readback](qfo_order_replay_readback_22328.json)
binds the scheduler, exact commands, inputs, software, logs and output files.

| Comparison | Groups per arm | Shared groups | Genes in changed groups | Unequal constraint positions |
| --- | ---: | ---: | ---: | ---: |
| Historical vs retained order | 351,739 | 351,739 | 0 | 0 |
| Historical vs canonical order | 351,739 | 351,737 | 49 | 5 |
| Retained vs canonical order | 351,739 | 351,737 | 49 | 5 |

All partitions cover the entire universe exactly once. All three traces
contain 40,169 valid directed constraints. The semantic comparison retains
constraint order and direction while ignoring gene order within each side
and numerical diagnostic fields. First unequal constraint position: 921
(zero-based). Each differing partition has two groups absent from the other;
49 counts all genes in affected groups, not necessarily 49 moved genes.

Retained-order candidate partition and raw merge trace hashes exactly match
the historical files. Canonical ordering changes both partition membership
and semantic constraints. This isolates an ordering effect on cached inputs;
it does not measure an accuracy gain or loss, nor prove upstream or downstream
equivalence. Historical QfO scores remain unchanged and must not be attributed
to the installed canonical end-to-end pipeline without further validation.

Next: trace these membership/constraint differences and freeze the necessary
downstream phylogenetic comparison before execution. Reuse unaffected expensive
artifacts only after validating their full dependency identities. No tuning or
retry was performed. Resource measurements are shared-host descriptions, not
controlled efficiency evidence. The DGX remains deferred.

Validation: 44 focused tests pass, including full-partition rejection cases,
direction/order-sensitive constraint comparisons, and changed arm provenance.
Independent readback checks pinned artifacts before and after the comparison.
Other publication gates remain open.
