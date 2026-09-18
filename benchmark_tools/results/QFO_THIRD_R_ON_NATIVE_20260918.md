# Third Original-Release Reconciliation-On Native Admission

Cell `p1_c0_r1` completed as task 21671_2 with exit 0:0 in 1:14:31.
Independent native admission 21673_2 completed with exit 0:0 in 1:37.
The committed native admission report SHA-256 is
`29131d4d9358a037696acd75abc519bf02f2447d261fe33cf5278540239424c8`.

All 976,504 original-release genes are preserved across 390,980 candidate
families and 393,774 root HOGs, with 731 split source families and zero
cross-source merges. Native output contains 4,975,326 pairs; its SHA-256
is `6be75fea965cb79c713b6e9720b286b52244fe9a08387024cc245a988c0afe27`.
Root-HOG integrity is not a pairwise accuracy result.

## Independent Recorded-Event Check

Applied the unchanged deterministic 64-family audit (seed 20260923) to
the 26,109 eligible reconciled families. The sample contains 1,321 genes,
2,572 nodes and 16,806 reconstructed pairs. All 16,806 match the native
sample exactly; no group-boundary filtering is applied in this C-off cell.
Report `qfo_event_pairs_p1_c0_r1_20260918.json` SHA-256:
`d9ad0c79c61ea406b76fa12e10aabe8f06d2ec9f2be820f0883cc5bcab444390`.
All 21 verifier tests passed.

This validates sampled saved-node-event to pair conversion, conditional on
the recorded mapping-conflict labels. It does not independently infer
trees, verify those labels, test biological orthology, cover bypass
families or establish a statistical error bound for unsampled families.
Do not compare independently selected family IDs across cells as paired
biological families.

Pair-preparation task 21675_5 started after native admission. Its terminal
result and complete pair manifest must pass before scoring is submitted.
The final reconciliation cell p1_c1_r1 is now running as 21671_3. All of
these are original-release experiments, not corrected-input results.
No interim parameter selection or uncertainty contrast is authorized.

## Pair Conversion and Scoring Submission

Task 21675_5 subsequently completed with exit 0:0 in 2:58. It retained
4,959,440 of 4,975,326 native pairs, excluding 15,886 pairs outside the
original-release mapping. The committed pair manifest
`qfo_factorial_pairs_p1_c0_r1_20260918.json` has SHA-256
`7dbdb5d226cc0993509b935e62426ea71852cea8cf750d9555ad0a393a41df27`.

Submitted scoring job **21711**, index 5, from frozen executor
`25f328d994765369cfae0382a21c3e7fdb3b7dab`, binding that exact pair-manifest
checksum. Independent admission job **21712** uses executor
`9680ccced0e351fa62e0e76c1f393a232d04d00e` with dependency afterany:21711,
index 5, and the same checksum. Both use the existing frozen wrappers;
submission is not score admission. Complete eight-cell paired uncertainty
still waits for all endpoints and validations.
