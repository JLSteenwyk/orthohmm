# Candidate Score/Order Control

This post-hoc mechanistic diagnostic follows the
[residual trace](OB_CANDIDATE_RESIDUAL_TRACE_20260926.md). It is not an accuracy
optimization, new benchmark endpoint or timing comparison.

Hold fixed the exact 62,885-group profile-refined seed partition, sorted gene
indexing, species membership classes, frozen installed scientific source and
the private 0.11.0 diagnostic environment validated for job 22320. Run only
the unchanged satellite_v2 candidate expansion with its frozen parameters.
No search, profile building, clustering, phylogeny or scoring is repeated.

Execute these five arms once, sequentially, retaining all outcomes:

1. Historical hit order, historical normalized scores, no self hits.
2. Historical hit order, fresh normalized scores, no self hits.
3. Fresh hit order, historical normalized scores, no self hits.
4. Fresh hit order, fresh normalized scores, no self hits.
5. Fresh original arrays including self hits, as a separate control.

The four factorial arms must have identical unique directed nonself keys.
Score alignment is by gene-pair key, preserving each chosen order. Values
are already normalized; do not normalize again or round them. The self-hit
control must retain the original fresh arrays and report the excluded count.
Historical/fresh species membership classes must agree up to label renaming.

The primary diagnostic outputs are exact candidate partitions and merge
traces. Validate full gene coverage and candidate seed/merge consistency.
Compare each arm to all others and subsequently to retained historical and
fresh candidate files. Exact reproduction of anchor outputs is a check, not
a precondition to conceal a contrary outcome. If the self-hit control differs,
retain that observation rather than treating the factorial as a complete
historical/fresh contrast. No score-selected seed or parameter change is allowed.

The prepared plan is
`benchmarks/work/ob_candidate_order_scores_20260926/plan.json`, SHA-256
`d53ec873a11e6bfff40f3cb8c83985fd67336e1535c98d765951349c3235e2ee`.
It records input, source, helper, runtime and environment identities. The
previously admitted pickle may be loaded only after its hash is checked.
One scheduler attempt, two CPUs, 32 GiB and 30 minutes; BLAS/OpenMP threads
one, PYTHONHASHSEED=0. Stop on error, preserve outputs and do not retry or
resume implicitly. Shared-host durations are descriptive only.

Preparation initially failed before creating the plan because the dependency
readback did not contain the historical cache record. The preparer now obtains
that identity from the independently hash-pinned retained search audit and
rejects conflicting identities across the two reports. No native attempt
preceded this correction. Fifty-six focused tests pass, including score/order
alignment, self-hit separation, duplicate/missing-key rejection and candidate
content checks. Native completion and admission must be recorded separately.
