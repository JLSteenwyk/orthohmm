# Native Overhead Panel Submission

Submission returned array21838 and the SSH process exited before the
main-host clock observation2026-09-18T23:23:30Z. All transfers and the DGX
non-executing selection preflight completed before submission. Launcher
75b7723 and recipe/authorization9c13e3a were committed and pushed first.

Requested partition=spark,node=spark-7ff0,exclusive20CPU/96GiB,
array0-17%1, one-hour task limit, no requeue, eligibility delay60seconds.
Native timeout remains900seconds. The18-task order, input identities,
commands, paired statistic and numerical budgets are frozen in
DGX_NATIVE_FRONTIER_OVERHEAD_PROTOCOL_20260918.md and its pinned plan.

Recipe SHA-256:
`64a05b5201e78a9d8d46f302879a49bd5cc470bf7813eb5d16e4c88d79c51edc`.
Authorization SHA-256:
`77875e0273883454c25fcb697916aedc10f5d0d3081baa77ead997d6f4aa0b47`.

No DGX SSH/SCP/remote log reads are permitted after submission completion
until the local Slurm controller confirms every task is terminal. Poll only
the local controller meanwhile. Retain every failure and flag; do not inspect
partial native measurements to change the plan or selectively restart tasks.
No scientific timing admission follows from submission or native exit status.
