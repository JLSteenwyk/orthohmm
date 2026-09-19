# DGX Measurement Identity Failure Diagnostic

The complete terminal panel was scanned with
`benchmark_tools.diagnose_frontier_identity`. All 5,955 retained snapshots
pass individual frontier validation. Observed between-snapshot changes
occur only in the two failed tasks:

| Task | Observations | Identity change |
| --- | --- | --- |
| 5 | point_0270 to point_0271 | `/system.slice/sysstat-collect.service` appears, device/inode `[32, 1171912]` |
| 5 | point_0271 to point_0272 | The same service scope disappears |
| 12 | point_0000 to point_0610 | `cups.service` inode changes from 1004500 to 1177260; `cups-browsed.service` from 1004596 to 1177356, device 32 unchanged |

The boot ID and benchmark target scope are unchanged across these
transitions. Individual snapshots are internally consistent; the failed
invariant is identity equality across observations. This explains the
reported `Boot, target or frontier identity changed` exception without
assuming a reboot or an OrthoFinder inference failure. Both failed tasks
retain native exit zero but remain failed measurements.

No service was stopped, counter repaired, gate relaxed, or task rerun.
Scope appearance/replacement does not quantify service CPU use or causal
runtime interference. Boundary-only observations cannot locate the change
within the intervening command. Even unchanged endpoints cannot exclude
scopes born and removed between observations. Interval CPU flags in the
16 validated tasks are a separate unresolved issue. The 27 scientific
scaling runs remain unadmitted.

## Evidence and Reproduction

[Machine-readable diagnostic and all snapshot hashes](dgx_pressure_identity_21889_20260919.json.gz):
gzip SHA-256 `a68ed02e12c94be505699963f6197e880f37817a4a3969d7716eafa385f47df5`;
decompressed SHA-256 `18ce0907d403aa2f5967adfd5972185fd359078adcab656f8ca8ca3d845ed5e3`.
The [complete panel audit](DGX_PRESSURE_OVERHEAD_AUDIT_21889_20260919.md)
remains unchanged.

```sh
/home/bizon/anaconda3/bin/python -m benchmark_tools.diagnose_frontier_identity \
  --archive benchmarks/work/pressure_overhead_archive_21889 \
  --output benchmarks/work/pressure_frontier_identity_21889.json
/home/bizon/anaconda3/bin/python -m pytest -q \
  tests/unit/test_diagnose_frontier_identity.py
```

The diagnostic completed successfully. Five tests pass: stable identities,
add/remove/replace classification, boot changes, invalid snapshot rejection
and terminal gating before native reads. The CLI requires a fresh output.
