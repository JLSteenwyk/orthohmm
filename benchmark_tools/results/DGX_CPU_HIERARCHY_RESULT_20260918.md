# DGX CPU Hierarchy Control Results

The [prospective protocol](DGX_CPU_HIERARCHY_PROTOCOL_20260918.md), probe and
hash-checking launcher were committed/pushed as967fb44 before execution.
Job21816 completed0:0 in6seconds, zero restarts, exclusive spark-7ff0,
20CPUs allocated to the node,2CPUs per batch task and2GiB requested memory.
The initial submission omitted partition=spark and was rejected before a
job was created. Adding the explicit node-compatible partition admitted the
job; no executed trial was retried or overwritten.

| Trial | Host busy CPU seconds | Job-parent outer CPU seconds | Native sleeping-step CPU seconds | Batch-step CPU seconds | Load localization |
|---|---:|---:|---:|---:|---|
| Quiet | 0.010000 | 0.008365 | 0.007011 | 0.001353 | Descriptive only |
| Completed burst | 0.800000 | 0.788677 | 0.001637 | 0.787040 | Expected pattern observed |
| Sustained batch | 3.130000 | 3.131226 | 0.011583 | 3.119643 | Expected pattern observed |

The completed load consumed0.750135312process-CPU seconds. The longer
condition consisted of four sequential0.75CPU-second workers. Batch counters
also include launch/exit and observer work. Job-parent minus summed step
increments were0.000001,0,0seconds. The sustained host-minus-job difference
was negative(-0.001226seconds) and remains signed, not clamped or corrected.

All original source hashes and read scopes/brackets pass local replay.
These observations support availability of disjoint native/batch counters
and detection of a known completed batch load in these controls. They do
not calibrate delay/overhead for process-heavy native pipelines or explain
the previous interval3flags. Parent/child counter agreement in three trials
is not a general bound. Non-CPU interference remains unmeasured.

## Evidence

- [Raw report](dgx_cpu_hierarchy_controls_21816.json), SHA-256
  `fd22c8a1bfb4fdd9d949438955920eb570173601cfa321ad0ccfef5d07dd24e7`.
- [Scheduler record](dgx_cpu_hierarchy_controls_21816_scheduler.txt).
- Local archive:`benchmarks/work/dgx_hierarchy_probe_21816/`.
- Remote retained recipe/results:`/tmp/orthohmm-hierarchy-probe.Gw4hMMmu/`.

47focused tests pass, including exact replay of all three actual controls,
source identities, load location and scheduler properties. Frozen historical
observer sources were not modified. No timing threshold changed, no timing
panel was rerun, and scientific timing admission remains false. The next
step is complete native-command integration of these hierarchy reads with
explicit accounting and overhead controls, not a speed ranking.
