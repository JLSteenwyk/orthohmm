# OrthoBench Initial Clustering Depends on Runtime

## Controlled Contrast

The preceding [weight control](INSTALLED_OB_CLUSTERING_PROBE_20260926.md)
produced the same partition under historical/fresh weights and an installed
runtime repeat. This diagnostic reuses its exact historical-weight payload,
ordered gene/edge indices, CPM 0.1, seed 4, retained isolates and environment
variables. Two predeclared calls use the current Anaconda launcher followed by
a repeat, each in a fresh process with one attempt and a 1,800-second timeout.
There are no reference labels, search, profile, refinement or phylogeny calls.

The [recorded historical replay command](ob_native_replay_v2_batch_20260916.sh)
used `/home/bizon/anaconda3/bin/python`. Its **current** environment reports
NumPy 2.2.6, igraph 1.0.0 and leidenalg 0.11.0, whereas the installed runtime
reports NumPy 2.2.6, igraph 1.0.0 and leidenalg 0.12.0. The frozen scientific
source root is explicitly inserted under Python `-I`. Imported `externals.py`,
`helpers.py` and `leiden_worker.py` hashes must match the installed modules;
igraph and leidenalg native extensions are separately fingerprinted.

This is a **whole-runtime contrast**, not a claim that the shared environment
has been unchanged since the historical replay, and not an isolated test of
the leidenalg version alone.

## Findings

| Initial partition | Groups | Repeat |
|---|---:|---|
| Installed runtime, leidenalg 0.12.0 | 58,827 | Exact in preceding control |
| Current Anaconda runtime, leidenalg 0.11.0 | 58,772 | Exact in this diagnostic |

All outputs cover all 251,378 genes exactly once. Across runtimes, 57,740
groups match; 1,087 installed-only and 1,032 Anaconda-only groups involve
23,352 genes. Both newly executed calls exit 0 and agree with each other.
The tracked [receipt](ob_clustering_runtime_probe_20260926.json) contains full
commands, input/output hashes, measured runtime identities and all comparisons.
Its SHA-256 is
`0cd3e0f86ed1033af1cb65f4759fce38052d1d2a996b1ce5c7e7c3f9f6559fff`.

The runtime change is sufficient to change initial grouping in this controlled
case. Tiny search-score differences were not sufficient under the preceding
installed-runtime control. Neither result assigns the final 0.284504-point F1
decrease to a specific cause. The 23,352 changed genes here are not an estimate
of final accuracy error or the causal fraction of the 30,881 candidate-stage
gene changes. Library versions, native builds and runtime loading remain
coupled in this contrast.

Next isolate the dependency change in a private environment before considering
any runtime pin or a downstream replay. Do not modify the shared installation,
promote an apparently favorable score, or redefine the frozen baseline.

## Verification and Reproduction

All 77 focused tests pass across source-identity gates and the preceding
payload, graph, hit and production-accuracy tests. Plans and native outputs
remain under `benchmarks/work/ob_clustering_runtime_probe_20260926`.

```sh
/home/bizon/anaconda3/bin/python -m benchmark_tools.probe_ob_clustering_runtime \
  --repo . --output /tmp/ob_clustering_runtime_probe
```

The destination must not exist. The wrapper verifies prior evidence, checks
scientific-source equality before native execution, records each result and
stops on failure without retry. Descriptive shared-host durations are not
comparative timing evidence. One repeat does not establish universal
determinism. Historical benchmark tables and scientific defaults are unchanged.
