# Native Hierarchy Quiet Control Results

Protocol/launcherscc3f63d and27-file verified recipebb4e8ae were pushed
before array21820. [Submission observation](DGX_QUIET_SUBMISSION_21820.md)
records the completed SSH session at21:45:58UTC and requested60second delay.
Actual first start was21:47:14UTC. All three tasks completed0:0 without
restarts, with exclusive20CPUs/96GiB and unchanged native commands.

Only main-host scheduler polls and local tests occurred between submission
completion and terminal confirmation. No DGX SSH/SCP calls were made during
that interval. The last task ended21:48:06UTC; archives were copied only after
the local controller reported all tasks terminal. This operator policy does
not establish that unrelated services or kernel activity were absent.

## Retained Outcomes

| Configuration | Native validity | Observed intervals | Flagged intervals, zero-based |
|---|---|---:|---|
| High-sensitivity | 98groups | 4 | None |
| satellite_v2 | 98groups,1835native pair rows | 8 | 3 |
| OrthoFinder full | 99checkpoint groups,1834native pair rows | 10 | None |

All before/after runtime/system/input identities, native output checks and
counter replays passed. All whole-command screens pass; the satellite_v2
interval flag remains. Its native-step residual is0.370167CPU-seconds
(0.370034average cores). Batch-step usage is0.002749CPU-seconds and
host-minus-job-parent residual is0.367417CPU-seconds. Parent-minus-summed-step
usage is approximately one microsecond. Negative residuals elsewhere remain
signed. No threshold, label or inclusion rule was changed.

Operator SSH/transfers are therefore not a sufficient sole explanation for
the persistent satellite_v2 flag. This does not assign the residual to any
specific process or establish an accounting defect. The disappearance of
other flags in one ordered array is not a causal overhead estimate or a
speed comparison. No further selective repetition is authorized here.

## Evidence And Limits

- [Validated result](dgx_hierarchy_quiet_smokes_21820.json), SHA-256
  `29d0268a44d77c23e4822a2a062cc9d1aa4fffe88279885e82aab0a32835441d`.
- [Recipe manifest](dgx_hierarchy_quiet_recipe_v1_20260918.json), SHA-256
  `53a24cacbfefed66589d89ed9ae1bf4af701aa3e90453f219c2712611c8efb45`.
- Raw archive:`benchmarks/work/dgx_hierarchy_quiet_21820/`,1112files,
  4640867bytes at audit, including terminal scheduler and original-input records.

84focused tests pass, including recipe byte identity, exact equivalence of
quiet/native launchers apart from relocated paths/report labels, native
counter replay and retained adverse evidence. The auditor is unchanged apart
from those path/report substitutions. No original result was overwritten.

The next unresolved timing question is attribution/accounting of host-minus-
job activity during native work. Native overhead, non-CPU isolation and
scientific inclusion policy remain unproven. These small-fixture controls
do not admit the original27timings or establish publication readiness.
