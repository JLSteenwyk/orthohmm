# Constructor-Format Results

Job21327 completed0:0 in13:53, using frozen executore919834. All workers
used minimal imports and the same saved graph, original edge order, multiplicity,
vertex count and weights. No optimizer or accuracy scoring was invoked.

| Format | Repeat | Mismatches Before Weights | Mismatches After Weights |
| --- | ---: | ---: | ---: |
| NumPy int32 | 0 | 6 | 6 |
| Python integer pairs | 0 | 0 | 0 |
| NumPy int32 | 1 | 0 | 0 |
| Python integer pairs | 1 | 0 | 0 |
| NumPy int32 | 2 | 0 | 0 |
| Python integer pairs | 2 | 0 | 0 |

Independent admission checked the fixed six-worker plan, frozen executor/runtime,
input format in plan/parent/snapshot/native records, minimal import identity,
278 file records, stage-witness consistency and complete post-weight endpoint
hash reconstruction. Original constructor arrays were intact in every worker.
The mismatch is the same six-edge pattern observed in preceding diagnostics.

Snapshot qfo_constructor_formats_verified_20260916.json SHA256
bfb9f49e1eb32afb31ab898d0101cc785da1b54110c70d7ea06ad0259fb87631.

## Interpretation And Follow-Up

All three Python-pair observations matched; two of three NumPy observations also
matched. This bounded sample is not a general correctness guarantee or a causal
diagnosis: formats change allocations and conversion behavior. The pre-weight
hash is implied by complete bounded witnesses, not separately recorded in full.
No production defaults or frozen benchmark environments have changed.

A checked-worker adapter is prepared for a separate initial-graph replay. It
requires this complete admitted format panel, exact saved payload files, one
recorded Python-pair constructor conversion, and the existing full endpoint/weight
gate before the unchanged frozen optimizer. Any graph mismatch aborts without
retry or partition selection. The post-optimizer graph must match too. Tiny-graph
tests exercise the real optimizer on intact graphs and reject a corrupted graph
before it can run; they do not prove large-graph behavior.

Next run bounded fresh checked replays, retain every partition and compare
repeatability without accuracy selection. Only after native integrity and replay
provenance are established can QfO factorial/parameter work proceed. Matching
historical counts alone will not establish full historical equivalence.
