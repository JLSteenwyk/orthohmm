# High-CPM Constructor Observation

Job **22119 completed 0:0 in 2:16**, using frozen executor
`e990f1d013a1baae8d753d1eb0a344aaa88f2d14` and the
[prespecified diagnostic](QFO_CPM_HIGH_CONSTRUCTOR_PROTOCOL_20260923.md).
The [retained report](qfo_cpm_high_constructor_22119.json) records one fresh
minimal-import Python-pair construction of the failed high-CPM stage's graph.

All 984,137 vertices and 25,501,180 edges were retained. There were zero
endpoint differences before and after assigning weights, including native
versus saved, constructor versus saved, and native versus constructor
comparisons. Native and saved ordered endpoint/weight fingerprints match:

- Endpoints: `f182b9e9b23f9158e1c36525bf546b58a4a0cf1f67928d3b829ab822833a0382`
- Weights: `6c978e4b078c167c33eb9c69fbedef4e239475e2d93c517a296a964806bcb042`

Post-completion review rechecked 261 retained source/input/runtime/observation
records and the result gate. The worker log is empty. GNU time reports
2:16.66 wall time, 128.34 user seconds, 8.28 system seconds and 3,176,380 KiB
maximum RSS. These are diagnostic shared-host observations, not comparative
timing measurements. No optimizer was imported or called, and no predictions
or accuracy scores were generated.

This establishes one successful construction, not the cause of the original
SIGSEGV. It does not reproduce the failed worker's complete import/allocation
history or exclude an intermittent failure. Both the original failure and
this success remain retained; no scientific output is replaced. Full replay
and partial-output reuse remain unauthorized. A subsequent diagnostic should
target the original frozen-import context with explicit phase markers rather
than treating this success as a reason for an automatic full-pipeline retry.
