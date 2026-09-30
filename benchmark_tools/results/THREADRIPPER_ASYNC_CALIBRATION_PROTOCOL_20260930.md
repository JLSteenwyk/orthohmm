# Prospective Asynchronous-Host Calibration

This is one new bounded engineering diagnostic on the approved Threadripper,
not a production scaling repeat. It is justified by the
[identified reporting failure in job 22378](THREADRIPPER_CALIBRATION_FAILURE_22378.md).
Original failed outputs, thresholds and source copies remain retained.
Do not stop other jobs/services or interpret a scheduler allocation as isolation.

The [prospective protocol](threadripper_observer_calibration_async_protocol_20260930.json)
has SHA-256 `7ed6720e659586f54c6f180aa69ecd08c51ec087c14593f03bef153e9a3a6b5d`.
It pins 827 harness Python sources, the unchanged private interpreter/controller,
the fresh submission script and a
[resource source amendment](threadripper_resource_endpoints_async_20260930.json),
SHA-256 `00a446ef9fecd2d54f8cd3efb45dee49294f1e233da8fbea6fcd38c231502cb7`.
The amendment binds the changed local collector and added observation-thread
helper; primary CPU/wall/peak-memory scopes and all 27 planned identities are
unchanged. No scientific inference setting or historical DGX helper changes.

Commit the protocol and tested source before one submission. Use a fresh output
at `benchmarks/work/threadripper_observer_calibration_async_20260930`.
One exclusive `bizon` allocation requests 64 CPUs/task, 128 GiB and 26 hours,
with no requeue. The workload remains 32 Python processes, four additional
threads each, 64 MiB retained/page-touched per process and 30 seconds of busy
work. Python threads share their process GIL; this is not 128-core scaling.

All original acceptance checks remain: at least five complete interior samples,
all witnessed thread identities/scopes/affinities, one-second maximum point
cost, 1.5-second maximum interior gap, contained worker CPU with the original
excess bound, and peak memory containing the retained 2 GiB allocations.
The existing report also retains its 0.5-1.5-second interval check. No gap is
discarded, no limit is widened and no overhead is subtracted. Host errors and
foreign CPU stay evidence, not approval for a quiet window.

The new observer thread starts before native release, takes periodic host
samples off the one-second native point loop, and is joined with a 15-second
bound before the final post-command host sample. Thread startup/failure/join
errors fail the attempt. This changes scheduling, not report schema or meaning.
The join/final scans remain inside the native step's 45-second parked release
window; real integration must verify that they finish in time.

After terminal scheduler state, inspect the native/report logs and run the
existing independent audit exactly once to a nonexistent output:

```sh
python -B -m benchmark_tools.audit_threadripper_observer \
  --directory benchmarks/work/threadripper_observer_calibration_async_20260930 \
  --protocol benchmark_tools/results/threadripper_observer_calibration_async_protocol_20260930.json \
  --protocol-sha256 7ed6720e659586f54c6f180aa69ecd08c51ec087c14593f03bef153e9a3a6b5d \
  --output benchmarks/work/threadripper_observer_calibration_async_20260930/independent_audit.json
```

Inspect numerical checks and accounting, not CLI exit alone. Failure remains
failure and requires separate diagnosis; no automatic replacement follows.
Even a pass would establish only sampled descendant/cadence/accounting checks,
not causal observer slowdown, full native environmental handoff, continuous
containment, host isolation, controlled timing or publication readiness.
173 focused tests and shell syntax pass before submission.
