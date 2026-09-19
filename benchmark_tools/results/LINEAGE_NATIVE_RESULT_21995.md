# Native Lineage Diagnostic Results

All corrected jobs 21995/21996/21997 completed with exit 0:0, followed by
recorder 21998. The recorder retained 407 polls, zero observation errors and
no missing jobs. Collection began only after all four jobs were terminal.
The earlier failed launches remain preserved; no selective reruns occurred.

## Archive And Validation

The complete transfer to `benchmarks/work/lineage_native_archive_21995/`
contains 109,645 regular files totaling 749,972,966 bytes: all native outputs,
the deployed recipe and the four-proteome input directory. The existing
auditor documented in `LINEAGE_NATIVE_AUDIT_READY_20260919.md` completed
successfully. All three tasks passed provenance, raw-measurement replay,
native-product validation and canonical output comparison with their pinned
same-method pressure-panel predecessors. No temporal-order or boot-domain
issues were reported. This does not establish biological accuracy.

Retained machine-readable audit:
`dgx_lineage_native_audit_21995_20260919.json.gz`, SHA-256
`0e3a77f785b645f44e0ca5cb8bbaeffbdb7106e2104350045e299c658de9050a`.
Complete controller capture: `lineage_native_scheduler_21995.tar.gz`, SHA-256
`5cd480f02dba6311955c42acaf36aeaf2814addd0adac68e0c5300de95b7d3e2`.
Large native outputs remain local, not committed. The audit includes their
file inventory and checksums; reproducing it requires these local archives
and the earlier pinned pressure-panel evidence.

## Diagnostic Measurements

These single-run measurements are descriptive, not admitted comparative
timings. Each method used the same 73,266 proteins and exclusive 20-CPU,
96-GiB allocation.

| Method | Native wall seconds | Native-step memory peak bytes | Observation points | Original flagged intervals | Narrow flagged intervals |
| --- | ---: | ---: | ---: | ---: | ---: |
| OrthoHMM high sensitivity | 549.565188408 | 2307461120 | 551 | 1 | 1 |
| OrthoHMM satellite_v2 | 808.475823379 | 2463068160 | 810 | 56 | 7 |
| OrthoFinder full | 612.115494517 | 7039598592 | 614 | 0 | 0 |

Memory is the measured native-step cgroup peak, not GNU-time maximum process
RSS. Both accounting forms remain in the audit with their distinct scopes.

All eight narrow flags are `excess_unassigned_cpu`. High-sensitivity interval
178 and satellite interval 178 have root-minus-job differences of 1,866 and
2,159 CPU microseconds, respectively, with narrow read overhangs of about
16.62 and 14.65 milliseconds. Satellite intervals 552, 699, 723, 728, 729 and
730 instead have differences of 192,286-300,567 CPU microseconds and read
overhangs of 1.04-3.93 milliseconds. These non-atomic signed differences are
not interference bounds or causal identifications. No flag is dismissed,
threshold changed, or wall time corrected on this basis.

## Remaining Requirements

The collector survived all three native workloads without the earlier
sibling-inventory failure. However, not all narrow intervals passed.
Environmental validity, scientific timing admission and publication readiness
remain false. The lifecycle control only tested completed user services;
during-read churn, the remaining CPU discrepancies, collector overhead and a
prospective inclusion policy still need resolution. The prepared boundary arm
has not been deployed. Do not promote these diagnostics to the final scaling
comparison or replace a failed complete panel with a successful subset.
