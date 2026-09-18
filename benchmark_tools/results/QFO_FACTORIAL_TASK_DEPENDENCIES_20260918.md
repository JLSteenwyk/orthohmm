# QfO Factorial Task-Level Dependencies

## Scheduling Change

Reconciliation task 21671_0 completed successfully while tasks 1-3 were still
running or pending. Previously, all native admission tasks waited for the
entire reconciliation array, and all R-on pair tasks waited for the entire
admission array. These barriers unnecessarily delayed scoring of completed
cells.

Updated the dependencies of the existing submitted tasks with `scontrol`,
without cancelling, resubmitting or duplicating jobs:

| Reconciliation | Native admission | Pair conversion | Scientific cell |
| --- | --- | --- | --- |
| 21671_0 | 21673_0 | 21675_1 | p0_c0_r1 |
| 21671_1 | 21673_1 | 21675_3 | p0_c1_r1 |
| 21671_2 | 21673_2 | 21675_5 | p1_c0_r1 |
| 21671_3 | 21673_3 | 21675_7 | p1_c1_r1 |

Each native admission now depends on `afterany:<matching reconciliation>`;
each conversion depends on `afterany:<matching admission>`. The existing
one-task-at-a-time throttle on each array is unchanged. `afterany` retains
failure inspection: native admission still rejects any unsuccessful upstream
task, and pair conversion repeats native admission before conversion. No
failed inference is admitted merely because a dependency became satisfied.

For each row the operational commands were:

```text
scontrol update JobId=<native admission> Dependency=afterany:<reconciliation>
scontrol update JobId=<pair conversion> Dependency=afterany:<native admission>
```

All eight updates returned success. A subsequent scheduler read showed the
expected task-specific dependencies and 21673_0 running on bizon. Its
materialized raw JobId is 21691; the public array task remains 21673_0.
Do not confuse scheduler materialization with a replacement scientific run.

## Scientific Invariants

No input, parameter, resource reservation, conversion semantics, executor,
reference, uncertainty protocol or result path changed. This workstation
correctness pipeline is not matched runtime evidence. Dedicated DGX timings
remain separate, and no unrelated jobs were stopped.

`admit_qfo_factorial_cell.require_success` requires the exact matching
reconciliation task to be COMPLETED, exit 0:0, on bizon with 32 allocated
CPUs. It binds raw scheduler IDs to execution provenance and checks successful
postflight evidence. `prepare_qfo_factorial_pairs.prepare` re-runs that gate
for R-on cells. These modules' focused suites pass 33 tests.

The checked, clean assessment executor remains
`25f328d994765369cfae0382a21c3e7fdb3b7dab`; the clean score-admission executor
remains `9680ccced0e351fa62e0e76c1f393a232d04d00e`. Scoring still requires a
successful terminal pair-conversion task and its explicitly frozen report
hash. No unscored output is counted as an accuracy result.

These are original-release experiments and retain the documented Xenopus
reference compatibility limitation. The corrected-release investigation and
separate rerun protocol are unchanged.

## First Task Progress

Native admission 21673_0 completed 0:0 on bizon in 1:38. Its report
`benchmarks/results/qfo_factorial_v1/native_admission_0.json` has SHA-256
`87aaf2236c7fbbff0aaa18f0d1f2e222b6d56a36d7402c630959b3172a3f024f`.
It verifies 4,966,346 native pairs, with no accuracy evaluation or score
admission. After excluding only the captured whole-array accounting string,
the report is structurally identical to the committed earlier independent
native review `qfo_factorial_native_p0_c0_r1_20260917.json`. No scientific
artifact changed between the two validations.

Pair conversion 21675_1 subsequently started on bizon under its existing
frozen converter. It remains subject to its own terminal completion and
pair-manifest checks before scoring. The full unit suite, including the new
corrected-input staging implementation, passes 2,518 tests in 52.31 seconds.
