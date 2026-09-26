# Reference-Size and Rare-Error Screen Results

Executed the six cells and 6,000 replicates in the
[protocol](SPARSE_DYADIC_F1_PROTOCOL_20260926.md), frozen in `e2a9eb0` before
execution. The [result](sparse_dyadic_f1_screen_20260926.json) has SHA-256
`c5231b5edec1aa052b8f73fc47811f3160b4ed149dc7706f2c58529570b4b405`.
The only real-data-derived simulation input was the reference-family
asserted-pair-count histogram; no tool-specific error rates were fitted.

| Case | Families | Coverage | Invalid variances / 1,000 | Mean estimated variance / empirical variance |
| --- | ---: | ---: | ---: | ---: |
| Unequal regular | 256 | 0.949 | 0 | 1.004 |
| Unequal regular | 16,844 | 0.947 | 0 | 0.956 |
| Unequal node | 256 | 0.933 | 0 | 0.939 |
| Unequal node | 16,844 | 0.941 | 0 | 0.926 |
| Rare perfect | 256 | 0.613 | 387 | 0.983 |
| Rare perfect | 16,844 | 0.651 | 348 | 1.040 |

The four regular cells pass the prespecified necessary screen. Both boundary
cells fail, so the broader screen fails. Invalid variances remain in the
coverage denominator; no clipping, continuity correction, retry, dropped
replicate or alternative interval was selected after inspecting results.
Eighteen related tests pass, including exact numerical comparisons between
sparse and full-grid variance calculations in each simulated regime.

## Why the Boundary Case Matters

With perfect recall in both simulated methods, their difference is determined
only by rare false positives. Both share a Poisson component; one has an
additional independent component of mean one. When that additional count is
zero, the observed paired contributions coincide and the plug-in variance
vanishes, although the population target difference is nonzero. That event
has probability exp(-1), about 36.8%, regardless of family count. The observed
invalid proportions, 38.7% and 34.8%, are consistent with this mechanism.
Increasing the number of families does not increase the fixed expected error
count in this stress design. Mean variance agreement alone therefore does
not establish Wald interval coverage or regular asymptotic behavior.

This does not prove failure for every native VGNC F1 contrast: the actual
methods have false negatives and are not these perfect-recall generators.
It does disqualify blanket use of the candidate in sparse near-boundary
settings. The earlier shared-clade dependence failure also remains; sparse
storage and reference-size matching do not repair it.

## Decision

No native confidence intervals are admitted. Further statistical work must
justify an explicit sampling target and dependence model and use an inference
method whose rare-event behavior is supported for the intended endpoint.
The current normal-Wald candidate is not promoted merely because selected
regular cases pass. Preserve VGNC point estimates, the eight-method count
decomposition and this negative validation in the manuscript. GO/EC, FAS and
the secondary mean remain separate uncertainty questions. This statistical
screen is not an evolutionary simulation or independent biological validation.

The protocol contains executable commands. The result records source/helper
hashes, NumPy version, seed, reference hash, all cell summaries and Wilson
Monte Carlo intervals. Existing result paths are refused. No benchmark score,
scientific default, installed tool or DGX state changed; the process completed.
