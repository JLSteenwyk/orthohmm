# Complete Native Hierarchy Results

Frozen launcher/auditor72fdb2c and verified transferred recipe1ee2db2 were
committed/pushed before array21817. All three sequential tasks completed0:0,
zero restarts, exclusive spark-7ff0 with20CPUs/96GiB. Scheduler elapsed times
were14,17,20seconds. All24transferred files matched local source bytes before
submission; recipe SHA-256 is
`f101970728f7e36340d5c7deb0682f2c1a0897d95b30bbe611a9ae7a66d61a1e`.

## Native And Counter Validation

All runs retained the frozen645-protein,eight-species inputs and commands.
Before/after runtime/system/input/enumeration identities matched, and native
output checks passed: high-sensitivity98groups; satellite_v2 98groups and
1835native pair rows; OrthoFinder99checkpoint groups and1834native pair rows.
These are validity counts, not benchmark accuracy or speed comparisons.

All complete-command hierarchy and original-threshold screens replay exactly.
High-sensitivity has4intervals with no flags; satellite_v2 has7intervals and
flags at3,5; OrthoFinder has10intervals and a flag at6 (zero-based indices).
All whole-command screens pass, without suppressing interval flags.

| Method/interval | Original host-minus-native CPU seconds | Batch-step CPU seconds | Host-minus-job-parent CPU seconds |
|---|---:|---:|---:|
| satellite_v2 / 3 | 0.338852 | 0.001652 | 0.337200 |
| satellite_v2 / 5 | 0.529818 | 0.001644 | 0.528175 |
| OrthoFinder / 6 | 0.644742 | 0.003046 | 0.641696 |

Job-parent minus summed step increments at those points are0,-0.000001,0
CPU-seconds. Thus unobserved batch-step CPU is not a sufficient explanation
of these recorded residuals. Host-minus-job differences still conflate
outside-job activity, asynchronous kernel work and accounting/read effects;
they are not rigorous foreign-load bounds or process attribution.

## External Observation Confound

The operator made an SSH log-read call during the array and copied the recipe
and input fixture from the DGX while the final task was active. Those commands
run outside the benchmark job and are a known potential competing workload.
Their individual CPU accounting was not captured, so they cannot be assigned
the residuals quantitatively. No quiet-host or overhead claim is valid here.
The integration remains useful for complete-command/native validation and
counter-scope checks, but these runs cannot establish scientific timing.

Any subsequent quiet control must be prespecified for all three methods,
perform transfers before submission or after terminal completion, and observe
the scheduler locally without SSH/SCP activity on the DGX during native work.
Retain this array and all flags; do not rerun only adverse methods or relax
thresholds. Such a quiet control is not yet executed or sufficient by itself
for scientific timing admission.

## Retained Evidence

- [Validated report](dgx_hierarchy_native_smokes_21817.json), SHA-256
  `51dc2656ca75d159363b913d863b6bd3eb093f5a6d4c3989c6a4f40585b0e735`.
- [Recipe inventory](dgx_hierarchy_native_recipe_v1_20260918.json).
- Raw archive:`benchmarks/work/dgx_hierarchy_native_21817/`,1108files and
  4604512bytes at audit, including terminal scheduler records and original inputs.

The first local audit began before SCP finished and failed on a missing
preparation file without writing a success report. After the same transfer
completed, the unmodified auditor passed. No native run was restarted and no
output was changed to satisfy the audit.79focused tests pass, including exact
retained-data replay, frozen recipe identity and adverse-flag retention.

No old timing is upgraded; overhead/accounting calibration, non-CPU isolation,
scientific inclusion policy and publication readiness remain incomplete.
