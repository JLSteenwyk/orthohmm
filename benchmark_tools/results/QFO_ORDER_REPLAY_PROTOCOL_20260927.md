# Paired QfO Candidate Ordering Protocol

Frozen before full execution on 2026-09-27. This is a reproducibility
experiment, not a new accuracy or controlled timing comparison.

## Inputs And Intervention

Use the admitted corrected QfO checkpoint (984,137 genes, 90,687,327 hits,
78 species), the same profile-refined seed partition, and the validated
installed recovery environment. Run `satellite_v2` once per arm, serially:

1. `retained_order`: original checkpoint hit order.
2. `canonical_order`: `canonical_hit_order_v1`, preserving hit values.

No search, profile refinement, phylogeny, scoring, tuning, or automatic retry.
The historical candidate output is a third, read-only comparison target.

Plan: `benchmarks/work/qfo_order_replay_20260927/plan.json`, 974,873 bytes,
SHA-256 `bbbb8b04e10d72839dc51de5b1b0484e2bc61a6802a189616cb64f09872c5c6b`.
The plan binds the launcher, canonical ordering policy, input and historical
output files, Python executable, and 3,052 installed wheel payload files.
Installed audit SHA-256:
`f1636e6e3222c135d2f174bb3f2c247f74e646d8653d087f068960f4926f2559`.

## Prespecified Readback

Require successful scheduler termination, both arm completion receipts,
unchanged pinned inputs/software, and verified output identities. Independently
validate full gene coverage with exactly one occurrence per partition.
Compare complete label-invariant partitions in all three pairwise contrasts;
report group counts, unmatched groups, and genes in changed groups.

Validate membership constraints against their candidate partition. Compare
the ordered sequence of source/target memberships, sorting genes inside each
source and target but retaining constraint order and direction. Numerical
diagnostic fields are not the semantic endpoint. Preserve raw outputs even
when their semantics match. Report every difference, including a mismatch
between the retained arm and the historical candidate output.

No differences establishes only candidate-stage equivalence on these cached
inputs, not end-to-end QfO reproducibility. Differences require investigation
before transferring historical scores to the installed canonical wrapper.
Neither result demonstrates an accuracy gain or justifies changing defaults.

## Execution And Validation

Request local node `bizon`, two CPUs, 128 GiB, four hours, no requeue.
Each arm has a 6,600-second timeout and an isolated interpreter with numerical
library threads fixed to one. Failed or partial results are retained; there
is no automatic retry. Shared-host runtime and GNU time maximum process RSS
are descriptive, not controlled timing or simultaneous process-tree memory.

Before freezing: 22 focused launcher/ordering tests passed. Both isolated
workers completed on a retained 16-gene, four-species fixture with three
groups and zero constraints. This checks plumbing, not nonempty merge-trace
equivalence. Fixture: `benchmarks/work/qfo_order_replay_fixture_20260927`.

Commit this protocol and launcher before submitting either full arm.
Record the actual scheduler ID and independent readback in subsequent receipts.
Other publication requirements remain open; the DGX stays deferred.
