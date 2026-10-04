# Generating-Tree Residual Mechanism Addendum

## Complete Error Cohort

The preceding [gene-tree oracle](SIMULATION_GENE_TREE_ORACLE_RESULTS_20261004.md)
and [upstream trace](SIMULATION_UPSTREAM_TRACE_REVIEW_20261004.md) retain all
70 cells / 10,125 candidates. Screen every candidate and select **all** with
within-candidate false positives or false negatives under the generating-root
arm: 20 candidates, eight tree-inferred eligible candidates and twelve bypasses.
This is a post hoc complete-cohort explanation, not an independently selected
test set. The [analysis plan](SIMULATION_ORACLE_RESIDUAL_PROTOCOL_20261004.md),
source `5a2cfcfc`, is pushed before new event-history tracing.

Trace all 800 cross-species pairs, including 380 true positives and 292 true
negatives. The 62 false negatives and 66 false positives exactly reproduce
the retained oracle's local counts. These selected-cohort counts do not replace
whole-panel accuracy or the much larger cross-candidate deficits. The
[generated tables](SIMULATION_ORACLE_RESIDUAL_RESULTS_20261004.md) are an exact
copy of the worker's Markdown; no candidate or negative result is dropped.

## Observed Mechanisms

All **62 false negatives** are true speciations that receive high-confidence
ortholog calls before native satellite-constraint filtering. Their genes lie
in the same generating-tree root group, but an unsupported satellite source
is detached and the pair is removed by the post-constraint membership filter.
They are not errors in the generating-tree LCA event call or root partition.
This isolates a specific recorded filtering decision on these fixed candidates;
it does not demonstrate that disabling constraints improves a real-data method.

All **66 false positives** have true duplication ancestors. Forty-five are
single-copy bypass predictions, where no native gene-tree inference occurs.
The other 21 receive speciation calls under the positive-paralogy rule because
the retained candidate has no overlapping species across that duplication.
For all 66, the full parent generating family **does** contain species-overlap
evidence at the same history ancestor. The candidate lacks that evidence.
Thus these observed cases differ from an original full family with complete
reciprocal loss: candidate retention removes an available overlap signal.
The analysis does not uniquely attribute that retention to search versus
clustering, expansion, or derived input masking. A single copy per retained
species is not proof of orthology.

These small residuals complement, rather than overturn, the preceding finding
that 99.765%/99.756% of divergent/divergent-turnover oracle false negatives
already cross candidate boundaries. Bigger or more sensitive candidates could
preserve useful duplication information, but indiscriminate merging could also
increase paralog contamination. No parameter or code default is promoted from
these development-exposed outcomes.

## Independent Check And Limits

The [compact readback](simulation_oracle_residual_readback_20261004.json) imports
no primary trace worker. It parses original reconciled XML, independently
induces candidate clades, maps them onto the species tree with Biopython,
recalculates fixed pair/root rules, reads native constraint events and applies
their support/filter semantics. All 800 pair ancestors, truth labels, stage
flags, overlap fields, classifications and exact original oracle counts match.
All 228 selected input/source identities rehash. Complete cohort selection and
summary arithmetic are checked against the pinned original oracle, not a list
chosen after tracing outcomes. This does not audit native HMM score calculations
or certify independent biological performance.

Initial readback `4594b045` fails before output because it resolves a
manifest-relative truth path instead of the recorded absolute path. Preserve
this execution. Corrected `55e7c9e3` adds four path-resolution tests; all 47
focused tests then pass and actual readback succeeds. The primary trace is not
rerun. Both successful diagnostic executions take seconds on the shared
Threadripper; these are not matched-resource comparative tool timings. No
unrelated workload is modified.

The [execution receipt](simulation_oracle_residual_execution_20261004.json)
pins sources, actual attempts, tests, local detailed output and selected public
summaries. The 1,462,878-byte detailed local report is not committed or described
as publicly deposited. Source, tests, plan, compact readback and generated tables
are committed/pushed at validated milestones. Prior timing, rc2, review PDF and
frozen benchmark scores remain unchanged.

This strengthens goal4.4 with event-history explanations. Full-goal uncertainty,
development-family inventory, label-independent biological strata, per-ablation
cost and publication-package integration gaps remain explicit. Next integrate
the validated gene-tree/upstream/residual addenda into the manuscript and
reproduction package while retaining incomplete scientific requirements.
