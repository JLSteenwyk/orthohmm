# Pressure Overhead Panel Released With Recorder

Deployment receiptfae1273 and source/protocol63d1098 were pushed before
submission. Submitted `run_dgx_pressure_overhead.sbatch` with `--hold` and
`--begin=now+2minutes`, receiving array21869. It requests18tasks, concurrency1,
exclusive spark-7ff0,20CPUs/96GiB, one-hour task limits and no requeue.
The native timeout remains900seconds. No scientific timing is admitted.

All DGX transfers, hash checks and non-executing preflight finished before
submission. Array submission was2026-09-19T04:23:38Z and eligibility was
04:25:38Z, leaving at least two minutes after preparation. Scheduler text
formats these times in local EDT, four hours behind UTC; retain the raw
record without rewriting its timezone presentation.

## Recorder Before Release

Exported the committed63d1098 controller capture source to
`benchmarks/work/pressure_overhead_capture_63d1098.py`, SHA-256
`c692cb0fc28084558fb4e4f2d21282a94a6a93ceb32255660bbb3e8d8d279c64`.
Submitted durable Slurm job21870 onbizon in the gpu partition,1CPU/128MiB,
25-hour scheduler limit, no requeue. Its24-hour capture deadline exceeds
the18one-hour task limits plus setup allowance. It polls the local controller
every5seconds; it does not read DGX files or launch remote commands.

The recorder initially waited for scheduler placement, then was confirmed
RUNNING and had saved its first successful poll before array release.
[First poll](dgx_pressure_overhead_first_poll_21869.json) SHA-256:
`e97b38bdbd8581daaa74520994bf7c6fd2f8119085a9a9adae7f1c74f31004f9`.
That record shows array21869 held, no inference allocation, and the requested
future start time. Released21869 only afterward, at approximately
2026-09-19T04:24:30Z. Local squeue subsequently reports BeginTime pending;
submission/release are not claims of native execution or success.

Recorder output is
`benchmarks/work/pressure_overhead_scheduler_21869/`; its job log is
`benchmarks/work/pressure_overhead_capture_21870.log`. Retain all poll files,
terminal `scheduler_INDEX.txt` records and final `capture.json`. A controller
capture failure remains a failure, not permission to synthesize records.

## Binding and Quiet Window

- Plan: `3950ccfa867c463ccfd7d8693dc331a3c85dc75ba67d9e8db89be4c1dda91c38`.
- Recipe: `11ce57a7895d92c456cc2b438157355cb0773cd473f4fbc4fc8305330a33eda3`.
- Authorization: `d1f5878580ed85829b1ce50ba0eebbcfd573975a6eb2b28471762a81978b9dea`.
- Submitted shell source: `b8f841d7c7fa90d792b4e37301d1df58addea2333b2223b5ae5be178e25d1544`.

[Array submission record](dgx_pressure_overhead_submission_21869.txt) and
[live recorder allocation](dgx_pressure_overhead_recorder_21870.txt) are
retained, not terminal-outcome substitutes.

No DGX SSH/SCP/log reads are permitted after release until local scheduler
state confirms all18tasks terminal. Poll only the controller during this
window. Do not use partial measurements to change task order, settings,
budgets or repetition counts. All historical21838failures and the original
27unadmitted scaling runs remain unchanged. Full analysis follows terminal
collection under the [frozen protocol](DGX_PRESSURE_OVERHEAD_PROTOCOL_20260919.md).
