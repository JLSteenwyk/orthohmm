# Complete Dual-Bracket Native Diagnostic

## Outcome

All three prospective jobs completed with scheduler exit 0:0, in order
21912, 21913, 21914. Recorder21915 also completed successfully. Collection
began only after all three native jobs were terminal; no DGX reads occurred
during their quiet window.

The complete archive audit validates all three command/runtime/input bindings,
all 437 archived recipe files, raw observation replay, native products and
canonical output equivalence with the prescribed earlier periodic tasks
1, 3 and 8. No clock-domain change or overlapping native command spans was
found. Full run inventories contain 109,169 regular files. All 73,266 input
proteins in four proteomes are accounted for. This is output consistency and
provenance validation, not an accuracy or controlled-resource admission.

| Diagnostic | Native command seconds | Intervals | Original flags | Narrow flags | Canonical outputs match |
| --- | ---: | ---: | ---: | ---: | --- |
| OrthoHMM high sensitivity | 550.250931349 | 551 | 327 | 0 | Yes |
| OrthoHMM satellite_v2 | 818.572113878 | 819 | 582 | 21 | Yes |
| OrthoFinder full | 621.868765903 | 622 | 177 | 1 | Yes |

These are descriptive engineering observations, not a speed ranking.
All original whole-command screens pass, but original interval failures
remain retained. All 22 narrow flags have reason `excess_unassigned_cpu`;
the largest flagged residual is 0.4178405458492093 average cores for
satellite_v2 and 0.2804104628400183 for OrthoFinder. No threshold, wall time,
failure status or output was altered. Narrow flags occur at satellite
indices 562,563,703,705,707,708,709,713,720,728,729,731,732,733,735,738,739,
740,741,742,743 and OrthoFinder index529 (zero-based).

## Native Products And Memory

- High sensitivity: 35,314 orthogroups.
- Satellite_v2: 35,627 orthogroups/root HOGs and 51,256 native pair rows.
- Full OrthoFinder: 24,052 checkpoint groups and 90,490 native pair rows.

Canonical partitions and applicable native pair sets match prior outputs.
The different product semantics must not be interpreted as equivalent
orthogroup counts across tools.

| Diagnostic | Native cgroup peak bytes | GNU-time maximum process RSS KiB | Native CPU PSI some microseconds |
| --- | ---: | ---: | ---: |
| High sensitivity | 2,266,378,240 | 569,272 | 1,944,399 |
| Satellite_v2 | 2,476,056,576 | 1,283,340 | 24,468,445 |
| Full OrthoFinder | 6,960,181,248 | 4,643,124 | 2,925,487 |

Cgroup memory includes wrappers/cache and is not maximum process RSS.
GNU-time maximum process RSS is not the simultaneous sum across processes.
PSI measures stalled time, not CPU-seconds or foreign interference. The full
report retains CPU full, I/O and memory pressure values, every interval,
both screens, GNU-time accounting and memory observations. Observation-window
endpoints are not exact native command boundaries.

## Evidence

- [Compressed complete audit](dual_native_audit_21912_20260919.json.gz):
  5,829,116 bytes, SHA-256
  `f6451d1fd5a96ff5f1a9b2155c6e4d9f17fabeb4de8a86cd53fbd1e70f37bff3`.
- Uncompressed work report: `benchmarks/work/dual_native_audit_21912_v1.json`,
  57,018,531 bytes, SHA-256
  `5139d9ac36b37962eab12064d872204fcbb69fb735ff562c542e297fbcf4c222`.
- [Terminal capture](dual_native_capture_21915_20260919.json): all three
  records retained, 438 polls, zero observation errors or missing jobs.
- [Independent poll replay](dual_native_scheduler_replay_21912_20260919.json):
  841 observations; no per-job gaps through first terminal detection, at
  polls107,294,437. Retained text and SHA-256 match those first observations.
- [Accounting](dual_native_accounting_21912_20260919.txt) preserves completed
  parent, batch and native-step records. Scheduler elapsed includes more than
  native inference and is not substituted for command wall time.
- Local archive: `benchmarks/work/dual_native_archive_21912/`, including
  frozen recipe, all native outputs, four original FASTAs and batch logs.
  Transfers contained 109,614 regular files and 1,029,003,788 logical bytes
  before adding local accounting/scheduler copies.

The audit ran with the command in `DUAL_NATIVE_ARCHIVE_AUDIT_20260919.md`
at revision `27c44b3`. Its source helpers are hashed in the report. Full
regression before collection passed 6,137 unit tests, native CLI integration
and all nine explicitly enabled legacy-runtime fixtures, documented separately.

## Interpretation And Next Gate

The full-command result does not reproduce an all-clean narrow CPU screen
for all three methods. Preserve every flag; do not launch a scientific timing
panel or quietly promote the narrow screen as an acceptance rule. Next is
mechanistic analysis of the flagged intervals using retained native-phase,
host/frontier and pressure evidence, with explicit uncertainty about attribution.
Any subsequent collector experiment requires a new prospective protocol.

Even clean flags would not establish overhead bounds, repeatability, larger
input behavior or non-CPU isolation. The repeated overhead experiment,
frozen scientific inclusion rules and 27 matched scaling runs remain open.
No historical timings, accuracy claims or publication status are upgraded.
