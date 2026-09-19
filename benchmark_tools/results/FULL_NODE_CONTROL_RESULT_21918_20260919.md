# Full-Node Control Results

## Execution And Replay

Job21918 completed0:0 in4:13 on spark-7ff0 with exclusive20CPU/96GiB.
Recorder21919 completed and retained47 contiguous controller observations,
zero errors and exactly one terminal observation at index46. Independent
replay reproduces that terminal text and SHA-256
`87265017944bc950396bd2336a934f2322da267e0aa268f3fffcfdcaf29e6ea0`.
Archive collection began only after both jobs were terminal.

Collected1,154 regular files,32,705,539 logical bytes, under
`benchmarks/work/full_node_control_archive_21918/`. All447 recipe files match
the pinned deployment. Before/after application/system/recipe runtime checks,
launch identity, allocation, prescribed trial order and embedded/raw evidence
validate. All nine trials independently replay both CPU screens and pass
worker identity, affinity, cgroup, duration, overlap, dose and exit checks.

| Trial | Condition | Original flags /21 | Narrow flags /21 | Common-work narrow flags /19 |
| --- | --- | ---: | ---: | ---: |
| 0 | Steady | 3 | 0 | 0 |
| 1 | Churn | 20 | 0 | 0 |
| 2 | Contended | 21 | 21 | 19 |
| 3 | Churn | 20 | 0 | 0 |
| 4 | Contended | 21 | 21 | 19 |
| 5 | Steady | 5 | 0 | 0 |
| 6 | Contended | 21 | 21 | 19 |
| 7 | Steady | 3 | 0 | 0 |
| 8 | Churn | 20 | 0 | 0 |

Each common-work subset uses complete enclosing original host windows, not
midpoints. Known outside-target CPU was detected in all three contended
replicates, including every common-work interval. All steady/churn narrow
intervals passed; original-screen flags remain retained, not erased.

Churn trials recorded192,482,198,466 and198,679 successful creation/wait
pairs, respectively. No creation cap was reached. Native worker self plus
waited-child CPU totaled397.492163,397.646906 and397.523554CPU-seconds.
Steady trials recorded398.896746,398.906472 and398.911804 worker CPU-seconds.
Competitor CPU was18.825855872,18.955834464 and18.953002432seconds, within
the prespecified5-21second range. These are actual process witnesses, not
inferred from requested CPU counts.

## Interpretation And Limits

The narrow screen separates these full-node native-only controls from the
known competing load. Process creation by itself, under this bounded
workload, did not reproduce the earlier21 satellite_v2 and one OrthoFinder
narrow flags. Those historical flags remain unexplained and unchanged.
This result does not establish absence of all interference, monitor overhead,
cause of historical flags, or eligibility of the27 scientific scaling runs.
No threshold changed and no scientific timing is admitted.

Prespecified residual distributions, native pressure, host context switches,
frontier accounting and within-block descriptions remain to be summarized
from the retained raw evidence. These correlated intervals are not independent
biological replicates and will not be used for interval-level significance.

## Artifacts And Tests

- Raw audit: `benchmarks/work/full_node_control_audit_21918.json`, SHA-256
  `a1a349d25d8a3758ac2e486d28b16a9d2a6febc1e10c2a15b243294065671f5e`.
- Compressed audit: `full_node_control_audit_21918_20260919.json.gz`, SHA-256
  `bd8e11ead6ee2e979784c8cd97730d3cf82ee62e4343c9cdf18c0475d2221cd3`.
- Terminal capture: `full_node_control_capture_21919_20260919.json`.
- Scheduler: `full_node_control_scheduler_21918_20260919.txt`.
- Audit entry point: `benchmark_tools.audit_full_node_controls`.

Focused suite123 tests passed, including nine new panel binding/scheduler
checks. Actual production replay above validates all nine real trial archives;
synthetic tests alone are not the evidence for this result. The earlier shell
launch failure21916 remains separately retained and is not a control outcome.
