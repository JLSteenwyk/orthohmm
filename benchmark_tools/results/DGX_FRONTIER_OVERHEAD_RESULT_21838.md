# Native Frontier Overhead Panel: Incomplete Evidence

Array 21838 is terminal for all 18 planned tasks: 15 completed and three
failed. Collection started only after the local scheduler terminal gate
passed. The complete archive is retained under
`benchmarks/work/dgx_frontier_overhead_21838/`; no remote measurement was
restarted or overwritten. This is an engineering overhead experiment, not
the 27-run scientific scaling panel.

## Independent Audit

The unchanged auditor completed successfully as a program, but did not admit
the panel. Its [full report](dgx_frontier_overhead_audit_21838.json) has SHA-256
`014837f097b9f907dba6f008c0eecc3b57a23e5365593444e52298c9a9db7d5d`.

| Outcome | Tasks | Count |
| --- | --- | ---: |
| Evidence validated | 0, 1, 2, 4, 7, 17 | 6 |
| Measurement wrapper failed | 3, 5, 6 | 3 |
| Required scheduler evidence missing | 8 through 16 | 9 |

Tasks 8 through 16 have retained `sacct` summaries instead of the detailed
`scontrol` records required by the frozen provenance audit. This was a
collection error. Completed records had expired from the controller when
recovery was attempted; summaries cannot authenticate all required fields.
The auditor's generic allocation-mismatch message is not evidence that those
tasks actually received different allocations. Original files remain intact.

Only high-sensitivity pair 0 is complete: output equivalence and the duration
gate pass, with periodic/boundary wall ratio minus one of
0.019860738842048035 (1.986%). This single pair cannot establish the planned
three-pair method median or complete-panel overhead budget. The other eight
pairs remain unavailable, not zero. No observed run is substituted for a
missing partner. All six validated arms pass their whole-command screen,
but this neither establishes isolation nor overrides unavailable arms.

The complete-panel numerical budget and environmental validity are not
established. Scientific timing admission remains false, including **0/27**
admitted scientific scaling runs. No overhead correction or controlled
method-speed ranking follows.

Reproduce the audit from the retained archive, at repository root:

```bash
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /home/bizon/anaconda3/bin/python benchmark_tools/audit_frontier_overhead.py \
  --archive benchmarks/work/dgx_frontier_overhead_21838 \
  --results benchmark_tools/results \
  --accounting benchmarks/work/dgx_frontier_overhead_21838/accounting_terminal.txt \
  --output /tmp/dgx_frontier_overhead_audit_21838.json
```

## Failure Diagnosis

The [wrapper failure summary](dgx_frontier_overhead_failures_21838.json)
preserves messages and hashes of the source verification and native-exit
records. Its SHA-256 is
`944e781f3631a59c3b20beb410f5c6cae4a470c754316e781d118f3cc5c20c01`.
All three native exit records report code 0 and no timeout; these are
measurement failures, not demonstrated native inference failures.

Examining every retained `point_*.json` frontier identity gives:

- Task 3: 799 points, one boot ID and one target throughout. Relative to
  point 0, only `/system.slice/apt-daily.service` is added, at points 6-36
  inclusive (31 samples), device/inode `[32, 1123383]`.
- Task 5: 614 points, one boot ID and one target throughout. Relative to
  point 0, only `/system.slice/NetworkManager-dispatcher.service` is added,
  at points 319-328 inclusive (10 samples), device/inode `[32, 1125007]`.
- Task 6 reports a frontier change inside a sampling operation. The failed
  operation does not retain both inventories, so its changing scope is not
  identified here.

For tasks 3 and 5 these additions explain the identity-check rejection in
`probe_cgroup_frontier.compare`; they do not quantify interference or prove
that the services changed native runtime. Neither service was stopped.

## Required Follow-up

Before another timing panel, make allocation-record collection automatic and
test it against controller retention expiry. Diagnose observation of transient
service cgroups without deleting their activity or weakening the current
panel's checks. Validate any revised measurement design separately and freeze
its inclusion rules before a new complete experiment. Do not selectively
rerun failed arms until the present panel appears to pass. Preserve this
negative result alongside any future panel.
