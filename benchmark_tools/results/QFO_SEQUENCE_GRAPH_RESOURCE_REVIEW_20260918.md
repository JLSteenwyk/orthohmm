# Corrected QfO Sequence Graph Resource Review

## Completed Prerequisites

Numeric source admission21792 completed0:0; both all-hit and top100
checkpoints match independently reconstructed source tuples. Payload
review21798 completed0:0 in00:05:01 with2CPUs/64GiB. Its retained
[report](qfo_graph_payload_20260918.json) has SHA-256
`d93b31c816fb1c05bfafe44265fe69d3276b64a5acd81acf1291b73e9e0c5e90`.

| Variant | Nonself finite directed hits | Named-array snapshot lower bytes | Snapshot upper bytes only | Input logical bytes, separate |
|---|---:|---:|---:|---:|
| All hits | 592530075 | 15735384968 | 22317930000 | 9500111012 |
| Top100 | 320184088 | 8926735267 | 13330512403 | 5142574804 |

These are verified dimensions and one source-derived array snapshot, not
peak-RAM predictions. Inspection of the frozen RBNH path confirms additional
winner/reciprocal/threshold/edge arrays; sorting, deduplication, checked
Python-pair graph construction, graph-library state, clustering and later
singleton/refinement work also require memory. The serialized parent/worker
workflow can retain more than one representation. None is captured by the
snapshot upper bound. Memory mapping does not guarantee resident savings.

## Prospective Allocation

Use32CPUs and384GiB for each variant, with the existing seven-day limit and
no requeue. This is a conservative allocation decision with substantial
headroom above the named snapshot and inputs, not a proven feasibility bound.
Use the same allocation for both arms; do not truncate hits to satisfy it.

At review, Slurm advertised192CPUs and1030000MiB onbizon, with64CPUs and
393216MiB allocated after the payload job completed. Linux reported
948339140KiB MemAvailable and4193276KiB SwapFree. One384GiB graph allocation
fits alongside the scheduled384GiB total; both graph variants must not run
concurrently. This availability snapshot is not a reservation or proof against
future load. Slurm scheduling and memory limits govern actual execution.

Run all hits first, then top100 with an afterany dependency. Top100 remains
an independently prespecified diagnostic even if all hits fails; do not
substitute its outcome, omit a failed arm, silently retry with new parameters,
or impute accuracy. Record failures and review any subsequent action explicitly.

## Frozen Execution And Admission

Generate the plan with the checked executor69adea4d5a1dc634783b20ff7bbe9fbe1c6464db
in`publication_qfo_sequence_graph_v1`. Its tracked benchmark/scientific code
was checked clean. Retain unchanged core7f3a9e4, BLOSUM62,CPM0.1,seed4,
profile expansion off, candidate expansion off and reconciliation off.
Require initial and multipass checked clustering calls and complete984137-gene,
78-species partitions at multipass and multipass_refined.

The plan's native command must not be run directly; use the checked driver
and batch entrypoint. Record the generated plan hash and submit only after
the plan/resource review is committed and pushed. Independent graph admission
must follow terminal execution and bind the actual result hash; conversion,
six-endpoint scoring and uncertainty remain separate gates. Neither allocation
nor job success alone admits scientific outputs.

This is incremental graph-only inference on the shared host, not dedicated
end-to-end timing or a matched-sensitivity claim. Existing HMM and Three
Kingdoms jobs are not stopped. No inference threshold or benchmark reference
is changed by this resource decision.
