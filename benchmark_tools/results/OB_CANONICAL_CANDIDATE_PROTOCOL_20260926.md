# Canonical Candidate Input Experiment

This experiment tests the benchmark-only `canonical_directed_pair_v1` policy
after the [order intervention](OB_CANDIDATE_ORDER_SCORE_RESULT_22322.md) and
[small regression fixture](CANDIDATE_ORDER_REGRESSION_20260926.md). It does not
change the frozen method, production defaults, score table or dependency pin.

Use the five retained factorial input definitions, fixed 62,885 seed groups,
gene indexing, species memberships, installed scientific source and private
0.11.0 runtime from job 22322. Immediately after constructing each arm's hit
arrays, sort unique directed hits by query then target. Keep scores unchanged
and retain self hits in the explicit full-fresh control. Do not renormalize,
round, alter thresholds, change seed groups or evaluate benchmark accuracy.

The adapter replaces only the diagnostic driver's input-construction function
inside the new process and restores it afterward. The retained driver and
all installed scientific functions remain unchanged on disk. The wrapper and
policy have separate source records in the plan; the native report continues
to identify the retained driver. This distinction must remain in readback.

Before inference, record dtype, shape and content hashes of each canonical
array and verify exact equality between historical/fresh order under each
fixed score vector. Verify that removing self hits from the canonical full
control gives the same arrays as the canonical fresh nonself arm. Different
score vectors must not be silently collapsed or rounded. Unit tests explicitly
retain this score contrast.

Run all five arms once, sequentially, preserving every candidate partition
and merge trace. Validate full coverage and reconstruct candidate seed/merge
consistency externally. Compare all ten arm pairs, plus retained historical
and 0.11.0 fresh candidate baselines. Retain unexpected differences; do not
select a policy using F1 or hide candidate changes from the existing defaults.
Full candidate agreement would establish only the tested ordering property,
not cross-version determinism, end-to-end agreement or improved accuracy.

The [preparation receipt](ob_canonical_candidate_preparation_20260926.json)
pins plan `benchmarks/work/ob_canonical_candidates_20260926/plan.json`, SHA-256
`6925ccc4684e1687a1f03ee8589d45d0f72b53605be5406f2a975c49acc8ecb7`.
A four-gene isolated fixture passed all five native arms and five external
readbacks. Forty-four focused tests pass. This is not full-data admission.

Allocate two CPUs, 32 GiB and 30 minutes; fix BLAS/OpenMP threads to one and
PYTHONHASHSEED=0. One attempt only, no automatic resume/retry. Timing is
shared-host descriptive evidence. Preserve job 22321's earlier startup
failure and job 22322's unmodified controls. Submission and terminal readback
must be recorded separately; no publication-readiness claim follows from
this diagnostic alone.
