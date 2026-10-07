# Accepted-Event Support Diagnostic

The [prospective protocol](NATIVE_QFO_CANDIDATE_SUPPORT_PROTOCOL_20261006.md)
was pushed in0a829b96 before selected feature distributions were inspected.
Exporter, independent rational-summary reader and focused tests were pushed
in48e7c391 before actual execution. Original inference, scoring, partitions,
aliases, uncertainty draws, defaults, manuscripts and rc5 remain unchanged.

The [actual report](native_qfo_candidate_support_20261006_v1/report.json),
[complete annotated pairs](native_qfo_candidate_support_20261006_v1/annotated_pairs.tsv)
and [unique implicated events](native_qfo_candidate_support_20261006_v1/implicated_events.tsv)
retain all40,690accepted events and all2,295changed VGNC pair paths. Each
accepted event contributes ONCE to descriptive summaries. The
[independent reader result](native_qfo_candidate_support_readback_20261006_v1.json)
reconstructs event/round/membership joins, every output row and rational
feature summaries, with exact identities/counts and1e-12absolute/relative
numeric tolerance. It imports neither exporter nor grouping kernels and
reuses, without replaying, the previous whole-partition/alias evidence.

## Descriptive Outcomes

| Event cohort | Events | Round 0 | Round 1 | Median support | Median source size | Median target size | Positive-infinity margins |
|---|---:|---:|---:|---:|---:|---:|---:|
| TP-only | 51 | 49 | 2 | 7.682283 | 1 | 17 | 12 |
| FP-only | 352 | 322 | 30 | 0.984776 | 9 | 32 | 55 |
| Mixed | 46 | 37 | 9 | 1.708556 | 1 | 43 | 10 |
| No changed scored VGNC pair | 40,241 | 33,699 | 6,542 | 0.273256 | 1 | 22 | 20,423 |

There are449directly implicated events, not2,183independent pair events.
Those events span150recovered asserted TPs and2,033added scored FPs.
The remaining12recovered TPs and100added FPs have transitive paths; all112
rows retain blank direct-event features. No direct event is invented for them.
All required serialized features are present. Positive-infinity margins retain
the original no-positive-alternative sentinel and are excluded from finite
min/median/mean/max summaries, not replaced with large finite values.

The TP-only and FP-only support ranges overlap:0.041871..24.914945 versus
0.030519..54.481907. Finite median margin is2.027995for TP-only and2.678003
for FP-only, respectively. These observations do not establish a discriminating
cutoff, significance, calibrated probability or a reason to change defaults.
Notably46events create both kinds of scored changes; their acceptable members
cannot be cleanly retained by simply labeling the entire event TP or FP.

This is development-exposed and conditional on ACCEPTED events, not rejected
alternatives or recomputed selection criteria. Reference scope, membership,
family size and native grouping confound cohort differences. The40,241events
without a changed scored VGNC pair are UNLABELED, not true negatives or correct
merges. A direct group-crossing event does not prove a direct pairwise HMM hit.
Neither feature association nor path tracing establishes a causal biological
mechanism, independent validation, valid VGNC confidence intervals, superiority
or publication readiness. No threshold fitting or parameter hill-climb follows
from this diagnostic.

## Execution And Anchors

Both actual commands ran once in isolated standard-library mode using the
retained scientific Python3.10.13. Export exited0at4.02s/280,432KiB maximum
RSS; independent reader exited0at21.03s/344,020KiB. Both recorded zero swaps.
These are shared-host postprocessing observations with unknown contention,
not native timing, repaired failure accounting or isolated-speed comparisons.
The [execution receipt](native_qfo_candidate_support_execution_20261006_v1.json)
retains actual arguments, observed resources and exact output anchors.

| Artifact | Bytes | SHA256 |
|---|---:|---|
| Report | 60,771 | a6dbe37a740ad657f9820cd20f78faf62c75e22e9a8ba28b258a886ee991ac8c |
| Annotated pairs | 843,587 | 07333b6ce3f7c563b08d084b9d7a514f0ab1c2d8d49885227284ac122ab35501 |
| Implicated events | 101,338 | 30d4c80d92e6d167d1b1fd7342302f46d50b55f982b6354bc287f8c708447c83 |
| Independent readback | 3,412 | dd399ccbf96943a67cb53061cafa2dd5f71cea3bc78538f4199a7770e2a4bf34 |

Original bound inputs retain alias report3b3ab8f1/readback4afc4aa6,
complete ledger3f0598f7 and original accepted traceec867165. Current bytes
are checked before/after each command. The old ledger is a report OUTPUT,
not a checked INPUT in the older report; its output anchor is verified directly.
No execution failure/retry, new native job, old evidence overwrite or rc5
rebuild occurred. Existing native23902and dependency review23910 remain the
same handles; their terminal outcomes and remaining cells are still pending.
