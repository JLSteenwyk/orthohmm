# Retained Host Pressure Diagnostic

Audited the host pressure-stall information (PSI) already retained in array
21838, without new native runs, changed admission rules or timing corrections.
The [machine-readable audit](dgx_pressure_21838.json) records 5,173 file hashes:
5,154 point snapshots, 18 native-exit records and the terminal accounting
inventory. It rechecks these hashes after analysis. The
[compact table](dgx_pressure_21838.tsv) includes every task and both memory
and I/O stall categories.

Audit SHA-256:
`2181a4a8aa01647e7aff55d5def46219570704905f118aaec1686a7903dd33f3`.
Source commit: `7e4a729`; source SHA-256:
`d17a31a22bcafeb69f8a267aab9775647ff1e6f50f5675552ceeab6d495af10b`.
Twenty new pressure tests and 34 existing panel-summary tests pass (54 total).

## Interpretation

The [kernel PSI documentation](https://docs.kernel.org/accounting/psi.html)
defines `some` as time with at least some tasks stalled, `full` as time with
all non-idle tasks stalled, and cumulative `total` in microseconds. System
CPU `full` is undefined and is not interpreted here. We difference cumulative
totals rather than treating the rolling 10/60/300-second averages as
measurements of the native command.

| Task | Scheduler state | Native window enclosed | CPU some (%) | Memory some (us) | I/O some (%) |
| --- | --- | --- | ---: | ---: | ---: |
| 0 | completed | yes | 0.137933 | 0 | 0.015079 |
| 1 | completed | yes | 0.323320 | 0 | 0.019453 |
| 2 | completed | yes | 2.900920 | 197 | 0.012141 |
| 3 | failed | yes | 3.104945 | 118 | 0.017434 |
| 4 | completed | yes | 0.518348 | 2 | 0.366117 |
| 5 | failed | yes | 0.570356 | 1 | 0.526013 |
| 6 | failed | **no** | 0.286787 | 0 | 0.025314 |
| 7 | completed | yes | 2.895463 | 227 | 0.011678 |
| 8 | completed | yes | 0.561068 | 2 | 0.504841 |
| 9 | completed | yes | 0.511636 | 1 | 0.347361 |
| 10 | completed | yes | 0.319475 | 0 | 0.020090 |
| 11 | completed | yes | 0.133211 | 0 | 0.014702 |
| 12 | completed | yes | 0.512978 | 1 | 0.348107 |
| 13 | completed | yes | 0.557683 | 2 | 0.572576 |
| 14 | completed | yes | 0.140077 | 0 | 0.011669 |
| 15 | completed | yes | 0.342701 | 0 | 0.024399 |
| 16 | completed | yes | 2.867333 | 173 | 0.013270 |
| 17 | completed | yes | 3.066845 | 177 | 0.015110 |

Percentages use the midpoint span of the outer host observations, not exact
native wall time or a share of CPU cores. Task 6 covers only 10.038107 seconds
before monitoring failed, not its entire native execution. Other tasks have
outer observations enclosing the recorded native window; that alone is not
continuous isolation evidence. Boundary-only arms retain two points.

All retained series have readable, nondecreasing resource totals and unchanged
host/observer identities. Memory `some` deltas range from 0 to 227 us over
these recorded windows. This is a description of host accounting, not proof
of zero memory contention or a bound on timing distortion. CPU and I/O stalls
include the native methods' own demand. Without native-scope PSI or independent
activity attribution, these data cannot distinguish external interference
from the workload being measured. No resource-based exclusion threshold is
chosen after examining these results.

The original audit still has six validated tasks, three wrapper failures,
nine missing detailed scheduler records and one available paired comparison.
Scientific timing admission remains false. Next measurement work must retain
native-scope and host observations separately, preserve transient diagnostics,
and prospectively validate overhead and inclusion rules; it cannot promote
this panel based on apparently small memory-pressure values.

Reproduce from repository root:

```bash
/home/bizon/anaconda3/bin/python benchmark_tools/audit_dgx_pressure.py \
  --archive benchmarks/work/dgx_frontier_overhead_21838 \
  --output /tmp/dgx_pressure_21838.json
```
