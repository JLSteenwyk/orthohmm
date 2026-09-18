# DGX Step Separation Positive Control

Protocol and implementation were committed/pushed as5461a20 before submission.
Job21800 completed0:0 in3seconds, zero restarts, spark partition/spark-7ff0,
exclusive20CPUs allocated,2CPUs/task requested,2GiB. Full terminal controller
text is retained in dgx_step_separation_scheduler_21800.txt; timestamps remain
unqualified. The raw trial report is dgx_step_separation_probe_21800.json.
Additional original files are retained in benchmarks/work/dgx_step_probe_21800.

| Fixed trial | Observer interval (s) | Host busy CPU (s) | Native cgroup CPU increment (s) | Burst process CPU (s) |
|---|---:|---:|---:|---:|
| Quiet | 1.000918807 | 0.02 | 0.007357 | Not applicable |
| Burst | 0.7861328545 | 0.80 | 0.005981 | 0.750070944 |

The observer/burst stayed in job21800/step_batch; the native waiting workers
were in step_0 and step_1. Native snapshots bracket both observer windows.
The burst exited successfully before the observer's final snapshot. The host
increment passed the frozen0.5 CPU-second positive-control threshold. All
snapshot fields were read without errors; retained raw replay, step separation
and monotonic brackets pass the added regression test. Scope labels in the
report have a doubled leading slash from joining path components; original
raw cgroup membership strings are retained and have a single leading slash.
No scope identity inference relies on the display label alone.

This is evidence that host totals retain a known exited batch child outside
a distinct sleeping native step. It is not an independent unrelated job,
general transient-detection calibration, overhead estimate or memory test.
Host totals include the observer and other work; no subtraction-based foreign
CPU estimate or statistical interval is justified by these two observations.
The native CPU deltas above include worker-side observation and waiting.

No scientific inference was repeated or unrelated process stopped. The
original27 timing runs remain descriptive-only. Next is a separately frozen
validation with compute/process-heavy native workloads and overhead controls
before choosing a new scientific observation/inclusion policy.
