# Threadripper Allocation Probe

Two bounded engineering jobs tested scheduler placement. Neither ran inference
or the resource collector; neither supplies comparative timing. Background
scientific workloads were untouched. The source was committed before each
job: `742f615` for 22338 and `2a6f6cd` for 22339. Nine focused tests pass,
including rejection of shared physical cores, wrong RAM caps, unlimited RAM,
incorrect allocation requests and changed child placement.

## Retained Failure

Job **22338 FAILED 64:0**, elapsed 00:00:00. A 32-CPU step with binding mask
`0xffffffff` was rejected before Python launched. Slurm reported allowed mask
`0x00000000000000000000FFFF00000000000000000000FFFF`, corresponding to
CPU IDs 0-15 and 96-111. This reserves both threads of 16 physical cores,
not the requested distinct 32-core placement. No probe output JSON exists;
the launch log remains `benchmarks/work/threadripper_allocation_probe_22338.log`.

## Successful Placement Diagnostic

Job **22339 COMPLETED 0:0**, elapsed 00:00:03, with an empty launch log.
The explicit corrective diagnostic reserved 64 scheduler CPU slots and used
the same 32-bit affinity mask. The [raw result](threadripper_allocation_probe_22339.json)
records parent and child affinity IDs 0-31, 32 unique package/core identities,
identical cgroup membership and the inherited allocation environment.
The effective RAM cap is **137438953472 bytes (128 GiB)** at job and step-user
ancestors. The child is a new Python process after a two-second sleep.

Both jobs reserved the node exclusively, so accounting reports AllocCPUS=192;
that does not mean 192 CPUs were used. For the successful diagnostic, distinguish
192 allocation slots, 64 step-reserved logical CPUs, and 32 affinity-allowed
native CPUs. The step-user cgroup cpuset is `0-31,96-127`; process affinity
is the narrower `0-31`. `cpu.max` is unlimited and `memory.swap.max` is `max`.
The RAM limit is not a combined RAM-plus-swap limit or a hard 32-CPU quota.

Exact successful submission from the repository root:

```bash
sbatch --parsable --job-name=orthohmm_threadripper_probe64 \
  --nodelist=bizon --partition=gpu --exclusive --nodes=1 --ntasks=1 \
  --cpus-per-task=64 --mem=128G --time=00:03:00 --no-requeue \
  --output=benchmarks/work/threadripper_allocation_probe_%j.log \
  --wrap='srun --exclusive --exact --nodes=1 --ntasks=1 --cpus-per-task=64 --cpu-bind=mask_cpu:0xffffffff /home/bizon/anaconda3/bin/python -B benchmark_tools/probe_threadripper_allocation.py --step-cpus 64 --output benchmarks/work/threadripper_allocation_probe64_20260928.json'
```

Do not rerun this command against its existing output path. This was a
deliberate engineering correction after terminal failure, not an automatic
benchmark retry or replacement of an unfavorable scientific measurement.

## Prospective Placement Clarification

For the [Threadripper amendment](THREADRIPPER_TIMING_AMENDMENT_20260928.md),
freeze native affinity to IDs 0-31 and retain 32 as each tool's explicit worker
setting. Reserve 64 scheduler step slots to accommodate this host's SMT mapping.
These identical rules apply to all three methods and all sizes/repeats. No
global topology, cgroup, service or scheduler setting was changed.

Before production launch, the local executor must distinguish requested worker
count from scheduler reservation, verify every observed native descendant's
affinity, retain effective CPU/RAM/swap limits, and fail admission on placement
drift. A process can explicitly widen its affinity within its allowed cgroup;
the simple inherited-child check does not prove every tool will preserve it.
Any monitoring gaps remain visible. Collector validation, runtime verification,
deterministic enumeration and a quiet window remain required. No timing run
is authorized merely by this successful placement diagnostic.
