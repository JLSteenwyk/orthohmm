# Matched-Graph Stage Trace Protocol

This is a post-hoc descriptive diagnostic, written after the final matched-graph
scores were inspected, but before stage-level outcomes were computed. It is
not a new confirmatory endpoint, parameter search or causal intervention.

Use all 35 reporting datasets and both search arms, exactly as admitted by the
strengthened graph readback. Do not rerun inference or choose representative
examples based on the results. Retain existing final-score definitions.

For each cell, enumerate unordered cross-species pairs at these stages:

1. Direct accepted search hits, collapsing either direction and excluding self
   and within-species hits. These are homology evidence, not ortholog calls.
2. Saved RBNH graph edges, and separately their connected-component closure.
3. Initial Leiden partition.
4. Saved singleton-augmented graph edges, and separately their component closure.
5. Multipass Leiden partition.
6. Final refined partition.

Compare every stage with the simulator's explicit ortholog pairs, not ancestral
homology-family pairs. Report TP/FP/FN, precision, recall and F1 as diagnostic
overlap statistics. For direct hits/edges/component closure, these statistics
do not imply that the underlying method predicts those pairs as orthologs.
Report gained/lost true and false pairs between adjacent partition stages;
report direct-hit-supported versus unsupported final true and false pairs.
Absence of a direct hit does not imply absence of an indirect graph path.

Aggregate metrics as equal dataset means within each condition and arm, then
equal condition means overall. Retain every cell and transition count. No
additional confidence intervals, significance tests or multiplicity-adjusted
advantage claims will be added. Verify final cell counts/metrics against the
already frozen score result before admitting the diagnostic.

Interpretation is stage localization only. Search identities, ranks and weights
all differ between arms; graph connectivity is not ground truth or a causal
explanation. Simulations remain development-exposed; five seeds share histories
across conditions. No change to the frozen method or publication defaults.
