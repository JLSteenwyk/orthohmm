# Retained Simulation Search-To-Candidate Trace

## Fixed Scope

Follow the completed gene-tree oracle result with an exploratory mechanism
trace under goal4.4. Prior outcomes are inspected/development-exposed. Freeze
the trace before calculating new stage counts: all70 generating-species-tree
OrthoHMM satellite_v2 cells, all retained reference ortholog pairs, seven
conditions and ten seeds20261101-20261110. Do not select bad seeds/families or
tune settings. Keep all pairs, including pairs successfully predicted.

Pin native admission17a1be71823b9fb7fb13082d05001d7667bb89c3366da2e3138a124e90907820,
tree preparation0caffab16f73019fabafe1fcfaac8c2de90960c5f1317803bf83308a6d2bf914,
and the preceding compact oracle readback
db016a1b4c68ab39caed061122e2196aadfb3ebdc38a4d4f5ad38d9edd664dec.
Check the selected reports, execution inventories, inputs/truth, numeric search
checkpoint, graph, candidate/seed/merge sidecars and native pairs against
recorded identities. Metrics must bind completed high-sensitivity built-in
search, the frozen parameters, and the selected native output directory.

## Observed Stages

For each true pair retain the presence of significant initial sequence-HMM
hits in both orientations, the direct final graph edge, undirected graph
component co-membership, candidate co-membership and native pair prediction.
The numeric checkpoint is written before graph inference and contains hits
whose native E-values passed the frozen cutoff, not every scored candidate.
Absence cannot separate prefilter rejection, scoring rejection, ranking/caps
or other search decisions. Scores are native sequence-HMM scores, not assumed
length-normalized values or new calibrated confidence.

The saved `orthohmm_edges.txt` is the final profile/singleton-augmented graph
written before final Leiden and refinement. It is not a saved intermediate
initial graph, native ortholog predictions or the final partition. Source
semantics are from frozen orthohmm.py SHA256
2afb89b9dc683e64e58208188f720e07701d4760ff09d1ac3a53c7a8075b84bb.
Use SciPy connected_components(directed=False), including isolates and
within-species edges. Paths can traverse genes outside the reference family.

Classify true pairs separated by candidates as:

1. Different retained graph components: no path in this saved graph.
2. Same graph component but different candidates: grouping/refinement/expansion
   leaves the pair separated despite connectivity.

Cross-tabulate each class by zero, one or two significant initial-hit
orientations and direct graph-edge presence. Also retain the full contingency
of all observed stages and counts for native false negatives within candidates.
Connectivity and direct hits are evidence, not interchangeable orthology
predictions or causal interventions. An absent direct hit may coexist with a
recovering indirect path; a direct hit does not require genes to share a group.

Validate candidate seed sidecar coverage, unique seed ownership and the frozen
merge-count arithmetic. Each candidate must be a union of distinct seed IDs.
The sidecar does not specify gene-to-seed membership inside merged candidates;
report this unavailable boundary rather than invent a pre-expansion partition.
No separate high-sensitivity run is silently treated as this native checkpoint.

## Reporting And Execution

Require exact truth-pair/native TP/FN reproduction and exact cross-candidate
counts against the preceding oracle readback. Count failures and unsupported
schemas explicitly; do not drop seeds. Keep a new all-pair TSV locally and
commit compact per-cell counts, source/input hashes and generated tables.
Report finite-panel counts/condition means without population intervals or
causal claims. No default promotion, native inference, new search or accuracy
benchmark rerun is authorized by this diagnostic. Preserve all prior outputs.

Use existing Python3.12.3/NumPy2.2.6/SciPy1.15.3/Biopython1.87, one numerical
thread, on the shared Threadripper without modifying unrelated work. Commit
the tested protocol/runner before executing the actual trace. Independently
read the TSV to verify all contingency/count projections and check graph
components with an existing independent library. This is not matched timing
or a new isolation/accounting gate.

The official [SciPy connected-components documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csgraph.connected_components.html)
supports the undirected sparse graph operation, not evolutionary truth or
causal interpretation. Remaining gene-tree oracle errors, genuine family
uncertainty, label-independent biological strata, development-family inventory
and per-ablation costs remain separate requirements.
