# High-CPM Frozen-Worker Boundary Diagnostic

## Prespecified Next Observation

The failed fourth worker of high-CPM job 22081_1 reached the Python-pair
adapter's before-constructor marker but not optimizer entry. Its saved graph
passed the exhaustive payload audit. Minimal-import and frozen-import direct
constructors (22119 and 22121) preserved all graph endpoints and weights, but
did not follow the complete frozen worker setup and constructor-adapter path.
Those successes neither diagnose nor rule out the historical crash.

Run one diagnostic subprocess on the exact saved graph through the frozen
worker, retaining the checked Python-pair adapter. Use a fresh payload with
read-only links to the four graph inputs and copied metadata with only its
output directory redirected to the diagnostic location. Preserve CPM=0.12,
seed=4 and include-isolates=true. Revalidate the failed-payload audit and
successful frozen-import diagnostic records before and after execution.

The worker must use the pinned launcher/import source, single-CPU affinity,
one-thread numerical settings and fatal-error tracing. Retain the worker's
module/library snapshot and constructor-adapter markers. Intercept the first
`leidenalg.find_partition` call without invoking the optimizer: compare every
graph endpoint against both saved arrays and the worker's actual `graph_edges`,
and compare the complete ordered weighted graph fingerprint with saved inputs.
Write and fsync an explicit preoptimizer-stop observation, then raise a unique
diagnostic termination exception. The parent must distinguish this expected
stop from any native signal, ordinary failure or unexpected normal return.

No optimizer, partition, groups, prediction or accuracy score may be produced.
Do not overwrite the failed payload, reuse partial scientific output, repeat
automatically, or proceed to a full native retry based on this diagnostic.

## Resources and Interpretation

Plan one CPU, 64 GiB, two hours, no requeue. This is shared-host diagnostic
accounting, not comparative timing. Success means this fresh execution reaches
optimizer entry with the exact graph. It is not the same allocation history
as the failed worker and does not identify an intermittent fault's cause.
Failure must retain traceback and last completed marker. Review either result
before changing the method or designing another experiment.

The `stop_before_leiden.py` observer is implemented with small native-igraph
tests, including endpoint/weight/constructor mismatch rejection and guaranteed
restoration of the patched optimizer function. Full parent/worker orchestration,
provenance binding and a frozen scheduler wrapper are still required. This
protocol does not claim that a full-scale diagnostic has run.
