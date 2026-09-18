# Frozen Initial Graph Trace

## Purpose And Scope

Reconstruct initial RBNH edges from the retained OrthoBench normalized-hit
cache using frozen core7f3a9e40dd7e79f842cc2c11fb8b548f9a802806. This is a
descriptive continuation of the all-family trace, not new search, clustering,
parameter tuning or independent biological validation.

The native function is observed at return to copy its endpoint thresholds;
it is not rewritten. Every reference pair's direct-hit scores are compared
against the admitted pair trace, and the threshold-derived decision must
match membership in the returned native edge set. Categories are no direct
hit, initial edge, no finite endpoint threshold, or below endpoint threshold.
Reference labels are read for pair reporting only after complete graph
construction. Cache insertion order, including its exact-tie behavior, is
preserved. The trusted pickle is checksum-verified before deserialization.

The accuracy.py source SHA-256 is
1a35944ab7fea859143f599b1262aa111272787f6f735b9acc37d299427f2ad6.
Admitted pair-trace report SHA-256 is
bda00fb593b357bc8f07e43544feae598150ae891ed82c988fb339a92349de0c.
Cache SHA-256 is
78a5af40ea2683a69e1baefefb0966c3549cffe71b24bda03b9ba4a6e29e65e1.
All are rechecked after reconstruction; the source and pair table are
recorded in the output report. Existing output directories are refused.

## Execution

Job21797 submitted from frozen executor
b84b69e5a66dcac34536d318cf4b1058e7e3d325 at
`benchmarks/work/publication_ob_initial_edges_v1`.
Controller confirms1CPU64GiB1h onbizon,no requeue,zero restarts; scheduler
subsequently reports RUNNING. Output destination:
`benchmarks/work/ob_initial_edge_trace_20260918`.
No result is admitted at submission. Terminal status, input/output identity
and pair/marginal arithmetic must be checked before manuscript interpretation.

29 focused tests pass:18 edge-trace tests and11 joint-count tests, including
eight randomized native graph fixtures with tied scores, exact-threshold
acceptance, asymmetric hits, infinite thresholds, invalid scores/thresholds,
edge disagreement and restoration of the tracing hook after failure.

## Limits

Initial edges exclude singleton assignment and later profile-added edges.
A missing hit cannot distinguish prefilter rejection from scoring rejection.
Presence or absence of an initial edge is not a causal explanation of final
membership because other graph paths and later operations remain possible.
Raw reference-pair memberships include within-species and low-certainty
pairs; they are not official weighted OrthoBench recall. Reconstructing this
graph does not independently establish historical runtime identity or provide
a matched-resource timing measurement.
