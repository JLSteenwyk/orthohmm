# Lineage Overhead CPU Flag Description

This post-outcome description uses the complete pinned
[18-task audit](LINEAGE_OVERHEAD_RESULT_21999.md). All nine periodic tasks
are included: 5,929 intervals, 23 narrow flags, 5,906 unflagged intervals.
The nine boundary tasks remain explicitly without interval coverage.
Original flags and all raw counters remain in the source audit.

Every narrow flag has reason `excess_unassigned_cpu`. The flagged intervals
have signed root-minus-system.slice differences of 186,298 to 1,228,942
microseconds (median 242,316). Corresponding system.slice-minus-Slurm-scope
differences range from -13,955 to 12,440 microseconds; Slurm-scope-minus-job
differences range from -11,124 to 900 microseconds. Negative signed values
are retained, not clamped or treated as invalid observations.

| Method | Periodic intervals | Narrow flags | Flagged root-minus-system median (CPU seconds) | Unflagged median (CPU seconds) |
| --- | ---: | ---: | ---: | ---: |
| OrthoHMM high sensitivity | 1,658 | 1 | 0.676781 | 0.001073 |
| OrthoHMM satellite_v2 | 2,416 | 20 | 0.241393 | 0.001651 |
| OrthoFinder full | 1,855 | 2 | 1.084799 | 0.000460 |

Flagged native average CPU use ranges from 0.200 to 19.492 cores; narrow
read overhang ranges from 1.042 to 5.595 milliseconds. These are descriptions
of differently timed measurements, not a decomposition or causal bound on
interference. In particular, the positive root-minus-system difference does
not identify a user process, root process, kernel work or delayed accounting.
The association with the narrow residual cannot establish its cause.

## Next Evidence Needed

The observation directs the next engineering control toward the root to
system.slice boundary. Prospectively observe the user.slice counter, root
direct membership and relevant host CPU categories alongside the existing
lineage reader under finite known workloads, with idle and native-only
controls. Separately test service creation/removal during a read, not only
services that finish between snapshots. Freeze workload ordering, limits and
evaluation before execution; do not stop unrelated services or workloads.

This is a next-step design requirement, not a submitted experiment or a new
timing-inclusion rule. No existing flags are removed, no thresholds changed,
no process is blamed and no scientific timings are admitted.

## Reproduction And Checks

```sh
python -B -m benchmark_tools.describe_lineage_overhead_flags \
  --audit benchmark_tools/results/dgx_lineage_overhead_audit_21999_20260919.json.gz \
  --output /tmp/lineage_overhead_flags_21999.json
```

The output path must be absent. The script requires the exact retained audit
hash and all 18 validated tasks, checks interval/flag agreement and lineage
path continuity, and retains signed complements. All 23 flagged rows and
all/flagged/unflagged distributions are in
`lineage_overhead_flags_21999_20260919.json`, SHA-256
`39cb0c7123887d15e54e54c099d9995042ba56a26c59e50cd4b13d0ec73f2bd2`.

All 50 focused description, lineage and distribution tests pass. Tests
cover signed values, boundary unknowns, unflagged coverage, missing/failed
tasks, contradictory flags, wrong paths/telescoping sums, nonfinite values
and wrong input hashes. This is an audited-data description, not a fresh
raw-measurement replay or a validated causal explanation.
