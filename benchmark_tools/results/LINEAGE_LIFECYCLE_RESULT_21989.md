# DGX Completed-Service Accounting Controls

## Execution

Three prospective controls completed in job 21989, exit 0:0, on spark-7ff0.
The exclusive allocation reserved 20 CPUs and requested 4 GiB; the observer
ran in a separate exact two-CPU Slurm step, with recorded affinity `[0, 1]`.
Each uniquely named owned service ran on CPU 0 outside the Slurm job and
exited before its final snapshots. No SSH polling occurred during the controls;
terminal status was checked from the controller before downloading evidence.
No unrelated service or existing scientific job was modified.

The [prospective protocol](LINEAGE_LIFECYCLE_CONTROL_PROTOCOL_20260919.md)
and deployed sources are frozen at `01b9c52`. All eight recorded Python-source
and protocol hashes were compared against that commit's Git blobs and match.
The transferred source archive was SHA-256 checked locally and remotely:
`df5ae398b77b4e3bfafff2f1be305feb233b71d269451cfb3ec6753d445808c7`.
Kernel: `6.11.0-1016-nvidia`, aarch64; Python 3.12.3. Executable hashes and
full scope identities are in [identity.json](lineage_lifecycle_21989/identity.json).

The preceding attempt 21988 failed its two-CPU affinity guard before any
control workload ran. An exclusive batch allocation exposed all 20 CPUs.
The corrected launch used a bound Slurm step and new source/output paths;
the workload and numerical response thresholds were not changed. Both
[accounting records](lineage_lifecycle_accounting_21989.txt) and batch logs
(`lineage_lifecycle_21988.log`, `lineage_lifecycle_21989.log`) are retained.

## Results

| Control | Burn process CPU (s) | User-manager aggregate delta (s) | Signed root-minus-observer delta (s) | Service absent before final reads |
| --- | ---: | ---: | ---: | --- |
| 0 | 0.750132 | 0.813623 | 0.813403 | Yes |
| 1 | 0.750131 | 0.814886 | 0.815919 | Yes |
| 2 | 0.750226 | 0.814222 | 0.817895 | Yes |

All three satisfy both fixed .5 CPU-second response checks. The service
membership was under the expected user-manager subtree, outside the observer
scope, with exactly the selected CPU affinity. In each trial the service
cgroup was already absent at the first post-exit existence check; that check
preceded the second aggregate snapshot. All raw before/after counters,
service stdout/stderr, removal observations, commands and trial reports are
retained under `lineage_lifecycle_21989/`.

The new raw-fixture test replays all six comparisons and all three assessment
results, verifies removal-before-final-read ordering, and independently
recomputes both response scalars from integer raw counters. Together with
collector and control tests, all 63 tests pass:

```sh
python -m pytest -q tests/unit/test_run_lineage_lifecycle_control.py \
  tests/unit/test_probe_cgroup_lineage.py tests/unit/test_probe_cgroup_frontier.py \
  tests/unit/test_probe_dual_cpu_brackets.py
```

Report SHA-256:
`aa6eef80de8666a83c711fb5b25b86aedf6e0cb2b9c067604cf40d7981cde5b7`.
Identity SHA-256:
`2e74a34dcb5f3b53250dbd3c6a5c4c0816a05666e60285e75f1294134efa430e`.

## Interpretation And Next Step

This supports the expected aggregate-accounting response after deletion of
these three owned user-service cgroups on this kernel. Other CPU activity and
service startup contribute to the deltas; they are not exact estimates of
burn CPU, causal interference or collector overhead. A pass does not establish
false-positive rates, behavior for every cgroup type, or absence of non-CPU
interference. The workload was deliberately introduced, not an uncontrolled
competing user job.

This control does not reproduce service churn during a counter read, the
specific system.slice services in the historical failure, or full-node native
workloads. Keep all historical failed points and residual flags unchanged.
Next integrate the reader as a distinct measurement schema with native
membership/host-bracket/pressure/memory checks and raw replay, then evaluate
the remaining lifecycle and overhead requirements under frozen protocols.
Scientific timing admission and the full publication goal remain unmet.
