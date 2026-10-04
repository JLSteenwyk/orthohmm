# OrthoFinder Four-Proteome Third Repeat Reviewed

This completes the [previously submitted attempt](THREADRIPPER_SHARED_SUBMISSION_22414.md),
index 18, full OrthoFinder 3.1.5, four proteomes, repeat 2. Fresh accounting
reports parent/batch COMPLETED 0:0 at 10:36 and native step COMPLETED 0:0
at 7:41. Waited for the parent post-runtime/environment checks; native-step
success alone was not treated as complete review. No inference was repeated.

The new bound terminal reviewer ran once and passed all four categories:
runtime, environment, resources and outputs-or-failure. Its
[public summary](threadripper_shared_attempt_22414.json) is byte-identical to
canonical review: 10,330 bytes, SHA-256
`889b0d10d7c298f2a127f8d7f989fa37db82d3aa7c7257dc46330ffab297b08a`.
Terminal controller is captured before purge. Runtime/input checks match the
repaired prospective binding before and after inference; raw resource and
whole-run process/pressure evidence are independently replayed. Process and
pressure failure maps are empty. Sampled foreign demand reaches 42.7759537979
CPU-core equivalents; preflight demand was 41.2694294498. Those observations
do not establish isolation or quantify causal slowdown.

| Primary Resource | Actual Value | Frozen Scope |
| --- | ---: | --- |
| Native wall seconds | 420.740000613 | Monotonic native-command interval |
| Native CPU seconds | 4455.373716 | Task-subtree CPU-stat bracket, including wrapper |
| Peak bytes | 7274196992 | Native-step lifetime, including launcher; not pure algorithm RSS |

Preparation, conversion, independent review and scoring are not added to the
primary native wall. Whole-job teardown CPU/peak endpoints remain unavailable,
not derived by subtracting overlapping scopes. The native validator checks all
73,266 input genes, 24,052 checkpoint groups and 88,890 expanded native pair
rows. These counts agree with the earlier two four-proteome repeats; counts
alone do not establish identical predictions, determinism or new accuracy.
No accuracy endpoint is evaluated here.

## First Complete Three-Repeat Cell

The unchanged new reporter checks actual nineteen-attempt history and direct
evidence, then produces
[snapshot v19](threadripper_shared_panel_snapshot_20261004_v19/panel.json):
262,270 bytes, SHA-256
`1204e49045fa5d8785dc34a61b3b01e57f75a3adcda39fe21b6adf3308e84c53`.
It retains 19 reviewed attempts, 18 with measurements, 17 eligible observations,
exclusions `[0, 17]` and eight not-yet-reviewed identities. Index 17 still has
null endpoints. All earlier eighteen rows remain unchanged.

Only full OrthoFinder/four proteomes now has three eligible repeats (22398,
22406 and 22414). Its wall median is 441.2889663s, range
420.740000613..474.299930051s; CPU median 4602.758499s, range
4455.373716..5008.896961s; peak median 7274196992 bytes, range
7268642816..7438495744 bytes. All other cell summaries remain unavailable.
These are descriptive observed ranges, not confidence intervals or a causal
efficiency comparison. No fastest-repeat selection or background subtraction.

The unchanged plotter produces
[figure v18](threadripper_shared_resource_figure_20261004_v18/shared_threadripper_resources.pdf)
with all 18 measured points, the retained gray-cross exclusion, the unplotted
abort explanation and exactly one complete-cell median/range per panel.
Actual PNG is visually inspected. All 52 focused reporting/continuation tests
pass in 60.14s, including PDF text bounds/panel pixels, all-three-repeat
arithmetic and missing-summary preservation. The initial added test incorrectly
assumed a standalone public 22398 receipt; correct it to use that attempt's
actual earliest retained snapshot. No result or production source changes.
Original failed test execution is not a failed scientific attempt.

Historical table/figure/archive/manuscript bytes remain unchanged. Supporting
raw/work evidence remains local, not a complete portable release. Full panel,
final manuscript/resource reconciliation and versioned/archive release remain
incomplete. Next frozen identity is index 19, high-sensitivity OrthoHMM,
four proteomes, repeat 2, only after full-prefix and fresh capacity/handoff checks.
