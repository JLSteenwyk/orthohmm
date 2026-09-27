# Candidate Ordering Regression Fixture

The [recorded fixture](candidate_order_regression_fixture_20260926.json)
exercises the retained candidate implementation (refinement source SHA-256
`991f1eb6a5f73d0442529ed19095a34b7c6ba8bff8dfe43a1127e24ec73fb26d`)
without altering production code or benchmark outputs. It follows the
[full-data order intervention](OB_CANDIDATE_ORDER_SCORE_RESULT_22322.md).

Five genes form a three-gene anchor and two singleton satellites. Twelve
directed hits supply bidirectional evidence with scores 0.1, 0.2 and 0.3.
Use one candidate-expansion iteration and cap attachments at one, with all
other arguments at the pinned implementation's defaults. This smaller cap
is a synthetic mechanism fixture, not the full satellite_v2 configuration.
There are no reference labels or claims of biologically correct membership.

Permute the three target-hit positions independently for each satellite,
applying each permutation to its reverse-direction hits as well: 6 x 6 = 36
orderings of exactly the same records. The unchanged function produces two
distinct partitions. NumPy 2.2.6's float64 reduction of `[0.1, 0.2, 0.3]`
gives 0.6, while the reverse gives 0.6000000000000001. Bounded candidate
attachment can therefore turn an order-dependent numerical difference into
a membership change in this fixture. This demonstrates a concrete mechanism;
it does not identify every arithmetic decision involved in the 93-gene
full-data difference.

The benchmark-only `canonical_directed_pair_v1` policy orders unique directed
hits by query then target under fixed lexically sorted gene indexing. It
preserves score values, dtypes and self hits, does not mutate inputs, and
does not normalize, round, threshold, deduplicate or change clustering logic.
All 36 fixture permutations produce the same candidate partition after this
canonicalization. Duplicate keys are rejected rather than silently given
a new interpretation. Empty hit arrays are allowed with a defined universe;
malformed indexing, nonfinite/nonpositive values and inconsistent shapes fail.

Forty-eight focused tests pass across the new policy, existing refinement
tests and factorial alignment tests. Source hashes confirm the worktree's
candidate implementation matches the frozen installation. The experimental
policy is not imported by production inference and is not a new default.

Next test the policy on all retained full-data score/order arms, including
input-array identity and candidate/merge readback. Preserve both historical
and fresh score vectors; do not select a policy using F1. If promoted as a
scientific change, use an explicitly versioned configuration and renew
end-to-end and independent validation. This fixture does not solve runtime
dependency sensitivity, cross-platform arithmetic, seed-group ordering,
score perturbations or the final phylogenetic F1 discrepancy.

```bash
python -m benchmark_tools.candidate_hit_order_policy \
  --fixture-output /tmp/candidate-order-fixture.json
python -m pytest -q tests/unit/test_candidate_hit_order_policy.py \
  tests/unit/test_refinement.py tests/unit/test_probe_ob_candidate_order_scores.py
```
