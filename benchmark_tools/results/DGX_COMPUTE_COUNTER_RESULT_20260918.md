# DGX Fixed-Work Counter Result

## Execution

Protocol and code were committed and pushed as `071989e` before submission.
Job21801 completed exit0:0 in52seconds with zero restarts on spark-7ff0,
spark partition, exclusive allocation, two CPUs/task and2GiB requested.
Slurm allocated all20CPUs because of exclusivity; each workload step requested
one CPU. This is not a20-CPU saturation experiment.

All six trials and24child processes completed. Each child performed the fixed
16MiB-buffer/256-hash workload and passed its checksum check. Native cgroup
CPU increments exceeded95% of summed child CPU in every trial; native memory
peaks exceeded the known buffer size. Every snapshot's optional counter reads
succeeded. The batch observer and native steps had distinct cgroup scopes.

| Trial | Monitoring | Work wall seconds | Child CPU seconds | Native cgroup CPU seconds | Native peak bytes | Observer CPU seconds |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| 0 | Off | 8.322864 | 8.133017 | 8.316948 | 34906112 | 0.002149 |
| 1 | On | 8.322638 | 8.132849 | 8.316231 | 34902016 | 0.021126 |
| 2 | On | 8.321992 | 8.133227 | 8.315907 | 34902016 | 0.021248 |
| 3 | Off | 8.320304 | 8.131299 | 8.314958 | 34906112 | 0.002166 |
| 4 | Off | 8.323710 | 8.133398 | 8.318259 | 34902016 | 0.002085 |
| 5 | On | 8.319261 | 8.130780 | 8.313553 | 34811904 | 0.020969 |

Monitoring-on trials each retained43batch snapshots, including the initial
snapshot; off trials retained only the initial snapshot. Both modes polled
the worker every0.2seconds. Median work wall times were8.322863894seconds off
and8.321992220seconds on, a descriptive change of-0.010473%. This is not
evidence of a speedup or proof of negligible monitoring overhead. Observer
CPU covers the polling interval, not all initialization/validation work.

## Evidence And Replay

Raw report: `dgx_compute_counter_probe_21801.json` (599588bytes).
Terminal controller evidence: `dgx_compute_scheduler_21801.txt`, preserved
verbatim including trailing whitespace. Complete per-trial handshakes/logs
are retained in `benchmarks/work/dgx_compute_probe_21801/`; the remote recipe
is `/tmp/orthohmm-compute-probe.eC4ElCVa/`.

The report records and replay tests check all three source hashes. The
source-pinned launcher also checked the frozen protocol hash before execution.
Current remote `/usr/bin/python3 --version` reports3.12.3. The local3.10.13
replay passed the same validation conditions; all summary fields matched
exactly except trial2's sum of child CPU times: retained8.133226752 versus
replayed8.133226751999999. The retained report was not altered. The regression
test permits an absolute1e-12 difference for this sum only, preserving exact
comparisons for all other fields. This is numerical replay, not bit-exact
cross-interpreter reproduction or a fully pinned system environment.

## Remaining Scope

This verifies coarse counter behavior with completed compute-heavy children
and measures observer cost on one fixed synthetic workload. It does not test
native orthology pipelines, multithreaded saturation, arbitrary memory peak
accuracy, thermal/frequency comparability, or unrelated-load detection.
The memory value is cgroup memory, not RSS. No host-minus-native estimate of
foreign CPU was admitted, and no scientific speed/memory rankings were added.

The previous27runs remain descriptive. A native pipeline smoke and a frozen
prospective inclusion/execution plan remain necessary before another scaling
panel. Both controlled_workload_verified and publication_ready remain false.
