# Installed OrthoBench Initial Graph: Score and Order Controls

## Design

Following the [retained-search comparison](INSTALLED_OB_SEARCH_COMPARISON_20260926.md),
cross historical/fresh normalized scores with historical/fresh directed-hit
order. All four arms use the same sorted gene indexing, species memberships,
frozen `build_rbnh_edges` implementation and audit runtime. Self-hits are
excluded; hit keys must match exactly and be unique. This is a post-hoc
mechanism diagnostic without reference-label scoring, not a new benchmark.
It reconstructs initial graphs only, without search, clustering or phylogeny.

Both production graph modules are checksum-verified. Inputs are pinned by the
preceding search receipt and rechecked after execution. The historical pickle
is admitted by checksum before deserialization. Numerical graph arrays and
thresholds remain in `benchmarks/work/installed_ob_graph_probe_20260926` for
subsequent diagnostics; their hashes are in the
[tracked report](installed_ob_graph_probe_20260926.json).
Report SHA-256:
`e930906950cc65a1b1141aa49b8e6b25d8b9d3e7ac15319505170bceb955c27f`.

## Results

Every arm has **1,803,122 undirected edges with exactly identical edge keys**.
At fixed scores, changing hit order changes neither edge weights nor gene
thresholds. At either fixed order, changing historical to fresh scores changes
1,031 edge weights and 120 gene thresholds, with maximum absolute difference
8.881784197001252e-16 in each. No threshold changes between finite and infinite.

| Controlled contrast | Edge additions/removals | Changed weights | Changed thresholds |
|---|---:|---:|---:|
| Order only, historical scores | 0 / 0 | 0 | 0 |
| Order only, fresh scores | 0 / 0 | 0 | 0 |
| Scores only, historical order | 0 / 0 | 1,031 | 120 |
| Scores only, fresh order | 0 / 0 | 1,031 | 120 |

All six pairwise comparisons, including the crossed arms, are retained in the
report. The unchanged topology rules out an initial edge-presence difference
under this common implementation/runtime. It does **not** prove that tiny
weight changes are irrelevant to Leiden or that historical native execution
used identical numerical libraries and gene/edge indexing.

The candidate partitions already differ before trees. The remaining interval
includes initial/multipass clustering, profile expansion, refinement and
candidate merging. Next test initial clustering with fixed indexing/seed and
both saved weight sets, then isolate runtime or indexing effects if needed.
Do not attribute the final 0.284504-percentage-point F1 decrease to roundoff,
Leiden, or reconciliation before a controlled downstream test.

## Verification

Sixty focused tests pass across the new alignment/edge-comparison tests,
retained-hit comparison, threshold tracing and production accuracy tests.
Input-array tests reject unequal/duplicate pair keys, self-pairs, nonfinite
scores and invalid endpoints. Comparison tests check reordered edges and
changed weights.

```sh
/home/bizon/anaconda3/bin/python -m benchmark_tools.probe_installed_ob_graph \
  --repo . --output /tmp/installed_ob_graph_probe
/home/bizon/anaconda3/bin/python -m pytest -q \
  tests/unit/test_accuracy.py tests/unit/test_trace_ob_initial_edges.py \
  tests/unit/test_probe_installed_ob_graph.py \
  tests/unit/test_compare_installed_ob_search.py
```

The output directory must not exist. No scientific defaults, historical
predictions or accuracy tables changed. This common-runtime diagnostic is
neither controlled timing nor a complete reconstruction of historical runtime.
