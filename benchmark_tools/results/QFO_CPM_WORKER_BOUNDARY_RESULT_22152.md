# Frozen-Worker Boundary Result

Job **22152 completed 0:0 in 1:46** on bizon, one CPU and 64 GiB, from clean
executor `d72ad55ab533763ac23e1fb86c511c67678e4a7b`. It followed the
[prespecified protocol](QFO_CPM_WORKER_BOUNDARY_PROTOCOL_20260923.md).
The [retained report](qfo_cpm_worker_boundary_result_22152.json) has SHA256
`376fc46a435f1988683630717c945f0f51751ab7f5fd5594178973bcad5e8be1`.

The frozen worker constructed the exact failed high-CPM graph through its
usual setup and checked Python-pair adapter. It reached the intercepted first
optimizer call with **984,137 vertices and 25,501,180 edges**. Native versus
saved, constructor versus saved, and native versus constructor endpoint
comparisons all found zero different edges. Ordered graph fingerprints match:

- Endpoints: `f182b9e9b23f9158e1c36525bf546b58a4a0cf1f67928d3b829ab822833a0382`
- Weights: `6c978e4b078c167c33eb9c69fbedef4e239475e2d93c517a296a964806bcb042`

The adapter records `constructor_returned` and the observer records
`stopped_before_optimizer`. The worker log is empty. Post-completion review
rehashed all 279 distinct path/size/SHA256 records from the result, including
source/input identities, worker snapshot, constructor marker and stop report.
The frozen executor revision and clean source were checked again.

GNU time reports 98.76 user seconds, 7.03 system seconds, 1:45.90 elapsed,
99% CPU and 3,185,076 KiB maximum RSS. These describe this diagnostic process
and children on a shared host, not controlled whole-pipeline efficiency.

## Interpretation

The graph and frozen worker setup/adapter path do not fail deterministically
in this fresh observation. This extends the earlier direct-constructor probes,
which did not use the full worker path. It does not reproduce all allocation
history, upstream validation activity or conditions of the historical SIGSEGV,
and does not establish its cause or exclude intermittent native corruption.

The optimizer was never called. No partition, groups, predictions or scores
were produced; high-CPM scientific outputs remain missing. Original failed
payloads and expensive prior stages remain unchanged. Any recovery must be a
separately reviewed checkpoint continuation retaining the failure and validating
reused upstream evidence, not silent replacement or an automatic full retry.
