# Complete-Command Interval Native Results

## Executed Integration

The [protocol](DGX_INTERVAL_NATIVE_PROTOCOL_20260918.md) and implementation
0808155 were committed/pushed before the first attempt. Array21807 failed
all three workers before native inference because the recipe omitted the
package initializer. The [packaging correction](DGX_INTERVAL_NATIVE_RECIPE_FIX_20260918.md)
was committed/pushed as f391e66 before a fresh array. Both attempts remain
archived; no failed output was overwritten or relabeled as a scientific run.

Array21810 then completed all three sequential exclusive20CPU/96GiB tasks
on spark-7ff0, exit0:0, zero restarts. Scheduler elapsed times were14,17 and
20seconds. Nineteen transferred recipe files matched committed bytes before
submission; the recipe manifest SHA-256 is
`d18b0c1dd2a81cdbe627d2aecb008005ad9a4de407529973a2be9dba6afe99cd`.
Before/after checks matched26,673 runtime and10,066 system inventory records,
the21-record recipe inventory, original inputs and native enumeration.

| Configuration | Groups | Native pair rows | Observation intervals | Flagged intervals |
|---|---:|---:|---:|---:|
| OrthoHMM high-sensitivity | 98 | Not evaluated | 4 | 0 |
| OrthoHMM satellite_v2 | 98;98 root groups | 1835 | 7 | 1 |
| OrthoFinder3.1.5 full | 99 checkpoint groups | 1834 | 10 | 1 |

All use the same645-protein,eight-species fixture. Native partition/pair/graph,
fresh OrthoFinder input-copy and GNU-time checks passed. These counts are
validity observations, not new accuracy scores. The command wall observations
3.090311123,6.755606178 and9.943522504seconds are engineering diagnostics;
do not derive speed rankings or overhead estimates by comparison to old smokes.

## Adverse Screens Retained

Complete native work was enclosed by observer endpoints and native before/after
reads. All intervals replayed exactly, and all whole-command screens passed.
However, the fourth interval flagged excess unassigned CPU for satellite_v2
(0.324569average cores) and OrthoFinder(0.269881average cores), above the frozen
0.25 threshold. Their whole-command residual averages were only0.058152 and
0.048555cores. No threshold was relaxed and neither run was selectively repeated.

These flags do not identify unrelated work: residuals can include kernel and
observer activity, cumulative-accounting delay and mismatched read windows.
Negative interval residuals were retained as well. Host outer windows overlap;
do not sum residuals or subtract them from native runtime. Explaining these
flags and calibrating accounting/observer effects remain necessary before
claiming controlled comparative timing.

Final native-step memory reads succeeded, with peak counters1005862912,
1007529984 and1089392640bytes respectively. All retained memory.events counters
were zero. These are cgroup peaks including wrapper/startup and cache, not
maximum process RSS or scientifically admitted memory rankings.

## Evidence And Tests

- [Validated result](dgx_interval_native_smokes_21810.json), SHA-256
  `4244bf36a04a01634cda3988eb5ed330ce7124400c2e8a347fc433ae1209659f`.
- [Failed attempt](dgx_interval_native_failed_21807.json), SHA-256
  `1a29b12c45baec472170172f8d917677c0fde52e118d28452e8b766365c0b8b4`.
- [Fresh recipe inventory](dgx_interval_native_recipe_v2_20260918.json).
- Raw successful archive:`benchmarks/work/dgx_interval_native_21810/`,
  1103files,4476817bytes at audit. Failed archive:
  `benchmarks/work/dgx_interval_native_failed_21807/`.

The first local audit lacked the copied original fixture directory and stopped
without writing a success report. Copying that directory from the DGX enabled
its existing input-hash checks and the complete native audit. No native output
or analysis result was changed to pass validation.

Ninety-three focused tests pass, including copied-package resolution from a
competing checkout, command failures/timeouts, invalid boundaries, exact raw
screen replay, failed-attempt retention and adverse interval retention. Output
validation uses the existing native validators, not a new independent oracle.

This completes short-fixture native-command integration, not general overhead,
accounting-error calibration, non-CPU isolation or a scientific run/repeat
policy. No27-run timing panel was relaunched; previous timings remain descriptive.
Scientific timing admission and publication readiness remain false.
