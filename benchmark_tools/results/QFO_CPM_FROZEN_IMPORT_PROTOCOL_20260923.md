# Frozen-Import Constructor Follow-Up

## Rationale and Fixed Scope

The [minimal-import diagnostic 22119](QFO_CPM_HIGH_CONSTRUCTOR_RESULT_22119.md)
successfully constructed the exact failed high-CPM graph, preserving all
25,501,180 weighted edges. It did not load the original OrthoHMM/Leiden modules
and cannot diagnose the historical SIGSEGV. This follow-up tests that specific
import-context difference, not another full scientific replay.

Run one fresh `frozen_imports` worker through the existing pinned
`probe_qfo_direct_graph.py` implementation. It imports `leiden_worker` from
`benchmarks/work/publication_qfo_replay_native_v1` and Leiden, verifies the
worker's resolved file, and replaces `leidenalg.find_partition` with a function
that raises if called. The worker directly constructs the undirected graph
using integer Python pairs and checks all endpoints before weight assignment
and all endpoints/weights afterward. No optimizer, partition or accuracy
analysis is requested. The frozen scientific source is not edited.

The parent additionally requires the completed minimal-import report with
SHA256 `40ca89f64fab8e7ecb88b069a2a45f26f9104dca812445a188d70b232032a25e`,
including its preserved source, runtime, inputs and observations. The original
failed payload is read-only; use a fresh output directory with only symlinked
input arrays. Python fatal-error tracing is enabled before native imports.
Retain log, exit code and last observation on failure; no automatic retry.

## Execution and Interpretation

Freeze the tested wrapper in a detached executor and submit
`qfo_cpm_high_frozen_constructor_20260923.sh`: 1 CPU, 64 GiB, 2 hours,
no requeue, one constructor only. Minimal mode remains the default; a
separate explicit mode selects this follow-up. Seventeen focused tests pass,
including native subprocess construction in both import modes and module
inventory checks. These small tests do not establish full-scale behavior.

Success would show that the frozen modules alone do not cause a reproducible
construction failure in this one observation. Failure would localize a
context-dependent symptom but not prove the imported module responsible.
The original worker's complete allocation/bookkeeping history is not replayed,
and neither outcome alone authorizes full replay, reuse of partial groups,
scientific score admission or a default change. Timing remains diagnostic
shared-host evidence; no DGX work or controlled efficiency claim is involved.
