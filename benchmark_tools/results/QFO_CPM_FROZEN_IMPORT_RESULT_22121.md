# Frozen-Import Constructor Result

Job **22121 completed 0:0 in 3:50**, from frozen executor
`0ca77ada17d5407de61205683aaf40782c2c3f02` under the
[prespecified protocol](QFO_CPM_FROZEN_IMPORT_PROTOCOL_20260923.md).
The [retained report](qfo_cpm_high_frozen_constructor_22121.json) has SHA256
`2e1aa3998e29c14c22babbd859dab45c5692d4f4e4aa1ec094982c3154c628c5`.

The fresh worker imported the frozen OrthoHMM/Leiden context, disabled
optimizer calls, and directly constructed the failed high-CPM graph using
Python integer pairs. All 984,137 vertices and 25,501,180 edges were retained.
Every endpoint comparison before and after weight assignment had zero
differences. Native and saved ordered fingerprints agree:

- Endpoints: `f182b9e9b23f9158e1c36525bf546b58a4a0cf1f67928d3b829ab822833a0382`
- Weights: `6c978e4b078c167c33eb9c69fbedef4e239475e2d93c517a296a964806bcb042`

Post-completion review recursively rehashed all 265 distinct path/size/SHA256
records in the report, rejecting conflicting identities for the same path.
The worker log is empty; the optimizer was not called. No groups, predictions,
accuracy scores or controlled timing evidence were produced.

This observation, together with the successful minimal-import diagnostic,
does not reproduce a deterministic construction failure merely from loading
the frozen modules. It does not identify the original SIGSEGV cause or
exclude an intermittent failure. The failed worker's full allocation and
bookkeeping history was not replayed. High-CPM scientific outputs remain
missing; no full retry or partial-output admission follows from this result.

## Low-CPM Continuation

Independent candidate validator **22086_0 completed 0:0 in 1:06**. Its report
is `benchmarks/work/qfo_cpm_candidates_admission_22086_0.json`, SHA256
`3a7602e1bde6fa4e8ba5f393700b112cbac9f63041e48b65578004c92af0a2bf`.
It admitted 984,137 genes in 361,165 candidate families from 397,045 seed
families and 35,880 recorded merges. This is candidate consistency evidence,
not accuracy or independent rescoring of search support. Recursively checking
both reports covered 3,303 distinct file identities without conflict.

The original whole-array dependencies would wait indefinitely for the failed
high arm. A [dependency-only amendment](cpm_low_dependency_amendment_20260923.json)
changed five pending low-arm tasks to successful same-arm predecessors:

```text
22088_0 afterok:22086_0  inferred phylogeny
22090_0 afterok:22088_0  native admission
22092_0 afterok:22090_0  pair conversion
22094_0 afterok:22092_0  scoring
22096_0 afterok:22094_0  score admission
```

Before/after snapshots verify unchanged array identity/throttle, CPU/memory,
time limits, no-requeue/restart settings, command and working directory.
No high-arm job, native command, scientific setting or admission check was
changed. Runtime gates still revalidate their own predecessor artifacts.
Phylogeny is resource-pending, not complete; failures stop this chain.
Shared-host resource contention remains unsuitable for comparative timing.
