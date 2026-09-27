# Installed OrthoBench: Initial Clustering Weight Control

## Scope

The [initial graph diagnostic](INSTALLED_OB_GRAPH_PROBE_20260926.md) found
identical edge keys with 1,031 tiny weight differences. Before native execution,
the diagnostic fixed three sequential calls: historical weights, fresh weights,
and one historical-weight repeat. All use the same 251,378 sorted gene IDs,
ordered endpoints, isolated installed runtime, CPM 0.1, seed 4 and retained
isolates. Each call has one attempt and a 1,800-second timeout; failures stop
the sequence without retries. No reference labels are read.

This uses the installed package's unmodified `orthohmm.leiden_worker`, not a
replacement clustering implementation. Runtime is NumPy 2.2.6, igraph 1.0.0,
leidenalg 0.12.0. The wrapper checks the frozen installation plan and source
records, records native extension hashes, limits numerical-library threads to
one, and launches each worker with Python `-I` and a clean environment.

## Result

All three calls exit 0 and produce **58,827 groups covering all 251,378 genes
exactly once**. Historical versus fresh weights yield no differing groups or
genes. The planned historical repeat also agrees exactly. All three output
files have the same SHA-256:
`20221c74398b2fc72cabcfcf4795512eb0ab51ab37f7dc90df87e4b8e9c609e7`.

The [tracked receipt](installed_ob_clustering_probe_20260926.json) includes
commands, execution records, payload hashes, runtime identities and strict
partition comparisons. Its SHA-256 is
`2fd69ad44dfe5ff5d1bb6f4d82b31523d6cce06d2a2e8173ad036ef68f527f13`.
Payloads, native logs, plan and partitions remain under
`benchmarks/work/installed_ob_clustering_probe_20260926`.

This rejects an initial-clustering effect of the observed tiny weight changes
**in this fixed installed runtime and configuration**. It is not proof that
roundoff never matters, that all runtimes give this partition, or that these
groups reproduce the historical intermediate partition. One fixed repeat is
not a general determinism study. It does not resolve the full-run score drop.

Next distinguish runtime effects from later multipass/profile/refinement and
candidate-merging effects. The historical graph construction/replay environment
is not automatically identical to the installed environment. Preserve the
same graph and seed for any runtime contrast; do not choose seeds or parameters
based on accuracy. No complete inference rerun is justified by this result
alone.

## Verification

All 74 focused tests pass across graph payload validation, directed-hit
comparison, graph comparison, native threshold tracing and production accuracy
helpers. After execution, input/source records and each partition/execution
record were rechecked before creating the compact tracked receipt.

```sh
/home/bizon/anaconda3/bin/python -m benchmark_tools.probe_installed_ob_clustering \
  --repo . --output /tmp/installed_ob_clustering_probe
/home/bizon/anaconda3/bin/python -m pytest -q \
  tests/unit/test_probe_installed_ob_clustering.py \
  tests/unit/test_probe_installed_ob_graph.py \
  tests/unit/test_compare_installed_ob_search.py \
  tests/unit/test_accuracy.py tests/unit/test_trace_ob_initial_edges.py
```

The destination must not exist. This executes the explicitly bounded three
native clustering calls, not a score-only reader. Wall times are descriptive
shared-host observations, not controlled comparative efficiency evidence.
No historical scores, defaults, main benchmark tables, search outputs, profile
expansion or phylogenetic predictions were changed.
