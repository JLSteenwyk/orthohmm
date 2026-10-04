# Native Upstream Trace: Mechanism Addendum

## What Was Tested

Following the [generating gene-tree diagnostic](SIMULATION_GENE_TREE_ORACLE_RESULTS_20261004.md),
trace all 163,527 reference ortholog pairs in the same 70 generating-species-tree
simulation cells. The [prospective protocol](SIMULATION_UPSTREAM_TRACE_PROTOCOL_20261004.md)
is committed with source `3d5fd9f7` before execution. All seven conditions and
ten seeds per condition remain, including successful predictions. No settings
are tuned and no native search, inference or timing measurement is repeated.

Use significant initial sequence-HMM hits, direct final graph edges, undirected
graph paths, candidates and native pair predictions. The saved graph precedes
refinement/candidate expansion; it is not native orthology or a saved initial
graph. Exact original native TP/FN/FP and the previous oracle's cross-candidate
counts reproduce. Seed sidecar/merge evidence is checked, but gene-to-seed
membership inside merged candidates is not invented.

## Results And Interpretation

The [machine-generated tables](SIMULATION_UPSTREAM_TRACE_RESULTS_20261004.md)
are an exact retained copy of the primary runner's generated Markdown.

In divergent simulations, 11,442 of 12,736 cross-candidate true-pair losses
(89.840%) lie in different retained graph components. In divergent-turnover
simulations the equivalent count is 11,554 of 13,072 (88.387%). All but one
pair in each of these component-separated categories lack significant initial
hits in both orientations. This is observed absence of retained significant
evidence, not proof of a k-mer prefilter failure, score-threshold failure or
absence of homology.

The other 1,294 and 1,518 losses (10.160% and 11.613%) are connected in the
saved graph but separated by candidates. They localize to the combined
grouping/refinement/expansion boundary. Respectively 117 and 128 even have a
direct graph edge. Paths may traverse other families, and an edge need not
force a group merge; this is not a causal estimate of Leiden or reconciliation
error. Candidate partitioning is fixed before the gene-tree intervention, so
generating trees cannot recover these cross-candidate pairs.

Outside the two divergent conditions, all cross-candidate losses are between
saved graph components and lack significant initial hits. Native false
negatives inside candidates remain in the tables, separate from these upstream
losses. These native within-candidate counts are not the smaller residual
within-candidate errors after the generating-root oracle intervention.

## Validation And Reproduction

The [independent readback](simulation_upstream_trace_readback_20261004.json),
source `72ad11eb`, imports no primary trace worker. It rereads every native
search-query/target array, graph, candidate partition, truth and native pair
table used by the TSV. All six stage flags, canonical pair IDs and ancestral
labels are checked against those actual inputs, and all 163,527 truth pairs
appear exactly once. igraph 1.0.0 independently reproduces the SciPy 1.15.3
component relations/counts. Every per-cell contingency, ancestral-family
count, total and condition projection matches. All 1,568 selected input/source
identities are rehashed. This is selected native data consistency, not complete
transitive runtime attestation or validation of the native HMM scores.

Forty-one focused tests pass. Both actual executions exit zero. The primary
trace takes 9.54 seconds and its independent readback 5.70 seconds on the shared
Threadripper. These are diagnostic process observations, not comparative tool
timing. No competing analysis is stopped or modified.

The [execution receipt](simulation_upstream_trace_execution_20261004.json)
pins commands, sources, outputs and tests. Detailed local report:
`benchmarks/results/simulation_upstream_trace_20261004/results.json`, SHA256
`eb096b6275b2035877bf67f7229c6983fd4f7010c25ba87620d558360fa78806`.
All-pair local TSV SHA256:
`012d9f395a6e3f6cdcc008caa9df7bdecd3b77ca9975b72f638ff0de15a9b62a`.
The compact readback/tables are committed; these larger local artifacts remain
pinned and are not silently described as publicly deposited. Prior rc2,
review PDF, original benchmark scores and defaults remain unchanged.

This strengthens goal4.4 stage localization but does not close independent
generalization, biological-error strata, family-level uncertainty, missing
development-family inventory, or per-ablation cost gaps. Next exhaustively
explain the remaining generating-tree oracle errors under original event
histories/constraints, then incorporate the validated addenda into a future
publication-package revision without rebuilding unchanged evidence merely
to create another milestone.
