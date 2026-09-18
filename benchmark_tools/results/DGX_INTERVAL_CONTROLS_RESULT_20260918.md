# DGX Interval Control Results

## Executed Controls

Implementation and [protocol](DGX_INTERVAL_CONTROLS_PROTOCOL_20260918.md)
were committed/pushed as f572260 before job21806. Slurm records COMPLETED,
exit0:0, zero restarts,21seconds on spark-7ff0. The exclusive allocation
reserved20CPUs; the batch task requested2CPUs and2GiB. Each native step used
one CPU. This was not an OrthoHMM/OrthoFinder benchmark run.

Both fixed-order controls met the frozen expectations. All22 observation
points were retained, yielding20 valid approximately one-second intervals.
No counter-read errors or invalid intervals were replaced or discarded.

| Control | Whole-window unassigned average cores | Maximum interval unassigned average cores | Flagged intervals |
|---|---:|---:|---:|
| quiet | 0.012356 | 0.060339 | 0/10 |
| completed sibling burst | 0.080175 | 0.779855 | 1/10 |

The sibling used0.750059216process-CPU seconds and completed between
observation points2 and3. Interval2, the third interval, flagged
`excess_unassigned_cpu`; its signed residual was0.780057CPU seconds.
The whole-window residual was0.801866CPU seconds over10.001387seconds and
passed the0.25-core average screen. Thus this actual control demonstrates
that a whole-window pass can coexist with a detected concentrated burst.
The difference between residual and injected CPU is not attributed solely
to foreign work: launch/exit, observer, kernel and accounting effects remain.

## Evidence And Verification

- [Raw report](dgx_interval_controls_21806.json),208702bytes,
  SHA-256`755b69b6837c54e0aada83aa581b7a1a0c39e3bd722270b5c875819f69a7e745`.
- [Scheduler record](dgx_interval_controls_21806_scheduler.txt),
  SHA-256`83040fa0a5bb78aebdd51c220c4574993d5caf0a5d754ab22600557770d6dc07`.
- Full copied recipe, per-point snapshots, worker logs and scheduler stdout:
  `benchmarks/work/dgx_interval_probe_21806/`. Original remote recipe:
  `/tmp/orthohmm-interval-probe.wtGflJd4/`.

The batch verified four source checksums and the protocol checksum before
execution. Sources were hashed again after execution. Local replay reproduced
every interval and whole-window result exactly from the raw snapshots.
Seventy focused tests pass, including synthetic invalid/missing evidence,
raw-result replay, burst timing/scope, source identity and terminal scheduler
state. The actual burst is specifically an excess-CPU flag, not merely any
kind of screening failure.

## Remaining Timing Gates

This validates one quiet and one burst control only. It does not estimate
false-positive/negative rates, bound accounting delay or show negligible
monitoring overhead across native pipelines. Outer host brackets overlap,
so interval residuals must not be summed. The fixed observation window is
not a complete native-command measurement; startup and teardown integration
is still required. Missing-read policy, non-CPU contention, memory scope,
thermal/frequency variation and scientific run/repeat/inclusion rules remain
to be resolved before a replacement timing panel. No existing timing value
was corrected or newly admitted. Scientific timing and publication-ready
flags remain false.
