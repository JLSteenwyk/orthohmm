# Native-Step Pressure Interface Verification

Added `probe_native_pressure.py` at source commit `401b20b`. It reads CPU,
memory and I/O pressure at the native Slurm step, not the observer's batch
step or one nested task. Raw per-resource counters and read timestamps are
enclosed by host observations. Validation checks job/step separation,
membership, cgroup device/inode identity, ordered read windows, raw/parsed
agreement and monotonic cumulative totals.

The pressure parser now accepts an explicit native-scope mode. System CPU
`full` remains uninterpreted by default; native cgroup CPU `full` is retained.
Historical host-pressure analysis semantics are unchanged. Its previous
source hash remains bound to the recorded historical commit, not the new
parser bytes.

## Verification

All 32 focused native/host pressure tests pass. Submitted local Slurm job
21865 on `bizon` with two CPUs, 256 MiB and a two-minute limit. The observer
ran in step_batch and the sleeping worker in step_0 with one CPU. Job,
batch and native step all completed 0:0 in the scheduler's one-second elapsed
resolution. This is not a comparative runtime or overhead measurement.

The [unchanged result](native_pressure_smoke_21865.json) records both raw
observations and five source hashes. SHA-256:
`ff0544e991f5124644cb1d5a2cee95c1db4e984ff6dd1f09cbde30f9d9a02cee`.
The [detailed scheduler record](native_pressure_smoke_21865_scheduler.txt)
has SHA-256 `e9cadc837e8f6fe28f0e247deb2a31b1a1dc7466b83823b32a7f89b5f5d1c4f6`.
A separate local replay reproduced every reported delta and verified all
five source hashes. Raw worker snapshots and logs remain in
`benchmarks/work/native_pressure_smoke_20260918/`.

Both observations identify native scope
`/system.slice/slurmstepd.scope/job_21865/step_0`, device/inode
`[29, 3156761]`, distinct from the recorded batch observer. Native CPU `some`
and `full` totals each increased 366 us; memory and I/O totals did not
increase. A sleeping worker still performs monitoring/handshake work, so
these are not expected-zero assertions or proof of an idle machine.

Command submitted from repository root:

```bash
sbatch --parsable --job-name=native_pressure_smoke --partition=gpu \
  --nodelist=bizon --cpus-per-task=2 --mem=256M --time=00:02:00 --no-requeue \
  --output=benchmarks/work/native_pressure_smoke_%j.log \
  --wrap='/home/bizon/anaconda3/bin/python -B benchmark_tools/probe_native_pressure.py --output benchmarks/work/native_pressure_smoke_20260918'
```

## Limits And Next Step

This verifies a local kernel/Slurm interface, not DGX portability, pressure
attribution or calibrated interference detection. The probe is not integrated
into a frozen inference benchmark. Native-step PSI includes descendants and
wrappers. Host and native pressure are not additive CPU-time counters and
must not be subtracted to infer foreign work. Non-atomic reads, accounting
delay and unobserved cache/bandwidth/thermal effects remain.

Before use in a scientific timing experiment, verify the interface on the
dedicated host, assess it with predefined workload controls, measure its
overhead, and freeze inclusion rules. Neither this smoke nor apparently low
pressure admits any earlier timing panel. No active inference executor or
historical output was changed.
