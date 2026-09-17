# DGX Collector Overhead Protocol

Freeze this panel before its paired outcomes. This is an engineering diagnostic,
not orthology inference or a tool performance comparison. No result authorizes
subtracting an estimated overhead from scientific timings.

## Calibration Evidence

The first fixed-per-thread workload21638 completed0:0 in4s, measured3.070083s,
averaging10.138466 cores over the observed span. Equal work leaves faster cores
idle while slower cores finish on this heterogeneous CPU. Preserve this valid
calibration; it does not establish sustained full load. Its collector report
SHA25627525ae6d1dad62e5a12d619fdeaddf1c6f5dbcc855c199cccf1f43ea671f803.

The revised fixture dynamically distributes deterministic10-million-iteration
chunks among20 threads. Each chunk owns a fixed seed/result index; scheduling
does not affect its checksum. Each thread uses a separate2MiB arena. Calibration
21639 completed0:0 in12s;1000000000 logical iterations per worker measured
11.887213s and19.580362 observed average cores. There were13 resource samples,
no errors, and maximum observed foreign CPU0.003351 cores. This is calibration,
not a paired collector comparison. Report SHA256
29add627cc1dd203ca01ce69ff9e33665dc0a0226b2e1ab20339591bc9127e03.
Both saved collector records passed independent replay of raw resources/host
observations. This does not certify absence of all external contention.

## Fixed Panel

Use spark-7ff0 exclusively through partition spark,20CPUs96GiB, CPU-only,
one task at a time,12-minute scheduler limit and600-second command timeout.
Run the same native workload with8000000000 logical iterations per worker,
20 workers, in each of six fresh allocations. This is approximately95s based
on calibration; the size is now fixed. No outcome-dependent resizing or retries.

| Array index | Pair | Mode | Resource wait interval | Host interval |
| --- | --- | --- | ---: | ---: |
| 0 | 0 | Sparse | 86400s | 30s, checked only after command-wait timeouts |
| 1 | 0 | Sampled | 1s | 30s |
| 2 | 1 | Sampled | 1s | 30s |
| 3 | 1 | Sparse | 86400s | 30s, checked only after command-wait timeouts |
| 4 | 2 | Sparse | 86400s | 30s, checked only after command-wait timeouts |
| 5 | 2 | Sampled | 1s | 30s |

Both modes retain the same wrapper, pre/post resource and host observations,
spawn/wait boundary, allocation checks and native output writing. Sparse is
not a completely uninstrumented run: it estimates the incremental cost of
periodic observation. Sampled uses the planned scientific cadence. Retain
all outputs and failures. Slurm array concurrency must be1, and exclusive
allocation does not prove that no unscheduled workload exists.

Source SHA2566da28892666f47053ae80ff66520897648f12a62413d739d918274c0f4c05cac;
ARM binary SHA25680f97405b9361fa021b3857298ee7c4a3ee69e2f99d774cfba55a37afb0851ec.
Compile command: gcc -O3 -fopenmp collector_load_fixture.c -o load, GCC13.3.
run_dgx_collector_load.sh SHA256
a323dacc0b0df4443b3dbdce3b1a786f6e4293ecac153de291e868d3a5aa8fca.
The runner checks fixture/binary/collector hashes before launch, uses the
existing isolated OrthoHMM Python environment, and sets OMP_NUM_THREADS=20,
OMP_DYNAMIC=FALSE, OMP_PROC_BIND=TRUE and OMP_PLACES=cores. Both versioned
remote recipes and calibration evidence remain retained.

## Admission and Reporting

Require terminal scheduler success, independent collector replay, exact
allocation/cpuset/memory cap, identical native checksums across all six runs,
expected20 workers and iteration count, unchanged recipe/runtime identities,
no observation errors, and at least60s per command. Sampled runs need at least
three successful host snapshots and mean CPU use>=18 cores. Report sparse
and sampled observation counts, observed host competition and sampling gaps.
Foreign CPU>=0.25 average cores in an observed interval flags contamination;
do not replace flagged runs automatically or discard their measurements.

Report each pair's sampled/sparse wall-time ratio minus1, its median and
range; also report CPU-time differences, observer scan durations and memory
measurements separately. No significance claim from three pairs. As a
predefined engineering budget, require median wall inflation<=5% and no
pair inflation>10%, in addition to the validation/load/competition checks,
before treating this bounded overhead diagnostic as acceptable. Negative
differences reflect runtime variability, not proof of negative overhead.
If a budget or gate fails, preserve evidence and revise the collector or
execution protocol prospectively; do not choose favorable repeats.

The fixture has one compute process,20 threads and roughly40MiB of worker
arenas. It does not cover large process inventories, disk-heavy workflows,
large-memory accounting, thermal stability on long runs, or every scientific
workload. Additional per-run observation costs and host warnings remain
visible. Passing this panel alone admits neither scientific outputs nor the
27-run scaling experiment.
