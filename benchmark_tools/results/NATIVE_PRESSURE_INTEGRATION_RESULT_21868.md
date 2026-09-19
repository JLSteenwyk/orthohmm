# DGX Complete-Command Pressure Integration

Exported422 committed Python files from8ffdc40 to a fresh DGX recipe. The
local and transferred tar archives match SHA-256
`d11bd9ddd20c82f8b14eac0a7532d991a742f45200cc74814d45c4e2968b1262`.
After completion, all422 remote file contents match archive member hashes.
The initial tar comparison reported UID/GID differences from extraction
ownership; a subsequent explicit byte-hash comparison verified contents.
No archive or source file was modified to suppress that observation.

Job21868 ran sequential periodic and boundary-only measurements, each with
`native_pressure=True`, around `/usr/bin/sleep 2`. The allocation requested
20CPUs/96GiB on spark-7ff0, five-minute limit and no requeue. Job and batch
completed0:0 in7scheduler seconds; native steps0 and1 each completed0:0 in3s.
Only controller polling occurred during execution; data collection followed
terminal accounting. No external service or job was changed.

## Measurements

| Mode | Points | Command wall seconds | Native CPU some/full microseconds | Memory/IO some/full microseconds |
| --- | ---: | ---: | ---: | ---: |
| Periodic | 4 | 2.042075876 | 2275 / 2275 | all0 |
| Boundary | 2 | 2.033524593 | 1978 / 1978 | all0 |

The scopes were native `step_0` and `step_1` under
`/system.slice/spark-7ff0_slurmstepd.scope/job_21868`, not the batch observer.
Their stable device/inode identities were[32,1150292] and[32,1150648].
Both reports bracket the entire native command. Raw per-point files match
the reports, and retained done records match native command results.
Independent parsing of raw PSI totals exactly reproduces all integer deltas
and finds no counter decrease. Both existing whole-command CPU screens pass;
this is not a validated isolation claim.

## Numerical Replay

Re-evaluation on DGX Python3.12.3 exactly reproduces both complete reports.
Local Python3.10.13 replay reproduces all pressure values, native measurements,
flags and other fields, but not three descriptive floating-point sums:

| Field | Local | Retained DGX |
| --- | ---: | ---: |
| Periodic interval0 root-minus-frontier CPU seconds | 0.0027380000000000043 | 0.0027379999999999974 |
| Boundary outside-target frontier CPU seconds | 0.006429999999999999 | 0.00643 |
| Boundary root-minus-frontier CPU seconds | 0.004559999999999995 | 0.0045600000000000016 |

The first strict local equality assertion failed on these differences.
They remain explicitly recorded as mismatches; no saved outcome, acceptance
threshold or replay implementation was changed to claim bit-identical
cross-runtime results. The differences are at floating-point roundoff scale.
The integration check establishes the pressure/command observation path,
not bit-identical cross-runtime summation of every descriptive field.

## Evidence and Reproduction

[Summary and raw-file hashes](native_pressure_integration_21868.json), SHA-256:
`58a38e75e74da7f8aa8744129bf8a75a4c183264b1e1267d32e631b7ae41a131`.
[Detailed scheduler record](native_pressure_integration_21868_scheduler.txt), SHA-256:
`c76ba44b9ce2150514419ba904e4b6a5f48e405ae722654f9ddd59d06e8312b9`.
Raw files reside in `benchmarks/work/native_pressure_integration_21868/`.
Source recipe: `/home/jlsteenwyk/projects/orthohmm-publication/pressure_integration_8ffdc40`.

For a new scheduled allocation with the same resource limits, import
`measure` from `benchmark_tools.measure_native_frontier_step` (periodic)
and `benchmark_tools.measure_frontier_boundary_step` (boundary), then call
each sequentially with fresh output paths:

```python
measure(['/usr/bin/sleep', '2'], fresh_directory,
        int(os.environ['SLURM_JOB_ID']), 20, 96 * 1024**3,
        60, 1., native_pressure=True)
```

For replay, call each module's `evaluate(report['points'], report['native'],
21868)` and compare with `report['screening']`; preserve the documented
local differences rather than silently rounding them away.

## Limits and Next Step

These are two short sleeping-command engineering checks, not inference
benchmarks or randomized repeated overhead arms. Do not calculate an
overhead ratio from their durations or infer absence of interference from
zero memory/I/O stalls. Pressure includes the native wrapper and descendants.
Scientific timing admission, controlled-workload verification and publication
readiness remain false. The next resource requirement is a separately
specified overhead experiment and defensible inclusion rules for actual
inference timings. Prior invalid panels and27unadmitted runs remain unchanged.
