# Leiden Distribution Isolates Initial-Clustering Runtime Difference

## Intervention

The [runtime contrast](OB_CLUSTERING_RUNTIME_PROBE_20260926.md) found different
initial partitions, but did not isolate a dependency. This diagnostic copies
the locally retained leidenalg 0.11.0 distribution into a new private overlay
and loads it with the otherwise unchanged installed clean Python environment.
No package manager operation or shared-environment modification is performed.

Eighteen files are copied, including package source, metadata, `_c_leiden`,
bundled libigraph and libleidenalg. Installed RECORD hashes are verified where
present; bytecode is excluded. Source/copy hashes are retained. Copy paths are
restricted to that distribution's three roots. Python runs with `-I -B`, and
only the private overlay is inserted before its installed package path.

Two calls (overlay and repeat) were fixed before execution, with one attempt
and a 1,800-second timeout each. Both use the same historical-weight initial
graph, ordered gene/edge indices, CPM 0.1, seed 4 and retained isolates as the
preceding installed control. The child records its imported module paths and
`/proc/self/maps` in the same process that invokes the unmodified scientific
worker. NumPy remains 2.2.6 and Python igraph remains 1.0.0. Scientific worker
and Python igraph native-extension identities match the clean installed run.

## Result

| Configuration | Initial groups | Agreement |
|---|---:|---|
| Clean runtime, leidenalg 0.12.0 distribution | 58,827 | Preceding control |
| Clean runtime, private 0.11.0 distribution | 58,772 | Exact current-Anaconda partition |
| Same private overlay, repeat | 58,772 | Exact first overlay partition |

All 251,378 genes are covered exactly once. The overlay recovers the complete
Anaconda partition, not just its group count. Relative to the unchanged clean
baseline, 57,740 groups match, 1,087 baseline-only and 1,032 overlay-only groups
involve 23,352 genes. Both calls exit 0. Readback confirms the child mapped the
three fingerprinted private native libraries, not the installed 0.12 extension.

This identifies the **leidenalg distribution, including bundled native
libraries**, as sufficient to reproduce the observed initial-clustering
difference in this controlled configuration. It does not identify an individual
upstream algorithm change or distinguish Python code from bundled native-code
effects. It does not establish the identity of the original historical runtime
or causally allocate the full-run F1 decrease.

The [complete receipt](ob_leiden_overlay_probe_20260926.json), SHA-256
`3cc5b1795d1d0b361a02e4d18fc98e0a1bc5307ab6c2c086fe28e7bc0e0c4e4a`,
retains copied-file identities, commands, module/native mappings, outputs and
partition comparisons. Local artifacts remain in
`benchmarks/work/ob_leiden_overlay_probe_20260926`. Native files are not added
to Git. The overlay is a diagnostic snapshot, not a recommended installation
or redistributable release bundle.

## Next Step and Verification

Use the same retained inputs and fixed configurations for a stage-resolved
downstream comparison before selecting a publication runtime pin. Do not
choose a distribution because it improves a development-exposed score. Keep
both results and distinguish reproduction of historical behavior from improved
orthology inference. No historical table, scientific default or dependency
installation was changed in this diagnostic.

All 88 focused tests pass, including overlay path restrictions and the previous
source, payload, graph, hit and production-accuracy tests. Post-run readback
checked native mappings against copied-file hashes and exact alternate-runtime
partition agreement. Existing source and copied-file hashes are rechecked.

```sh
/home/bizon/anaconda3/bin/python -m benchmark_tools.probe_ob_leiden_overlay \
  --repo . --output /tmp/ob_leiden_overlay_probe
```

The output must not exist. The command requires the fingerprinted local 0.11.0
distribution and the frozen clean environment; it does not fetch dependencies.
Shared-host durations are descriptive, not controlled performance measurements.
