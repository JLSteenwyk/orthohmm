# Native Frontier Engineering Result

Array21831 completed all three tasks with exit0:0 and zero restarts under
exclusive20CPU/96GiB allocations. Local scheduler start/end times (EDT):
high-sensitivity19:02:44-19:02:58; satellite19:02:58-19:03:15;
OrthoFinder19:03:15-19:03:35. No DGX remote calls occurred between submission
completion and local confirmation that all tasks were terminal.

The local replay audit validated frozen native commands, before/after
runtime identity, original/copied inputs, output partitions/native pairs,
raw observation agreement, final memory evidence, and both counter replays.
Archive: `benchmarks/work/dgx_frontier_native_21831/`,1117 files,8551104 bytes.
Result: `dgx_frontier_native_smokes_21831.json`, SHA-256
`2f6d33f9dad301bbe3d51d3135ca61ee883373e556b0cc852e7c66fa631c1fca`.

| Method | Native wall seconds | Native output | Intervals | Flagged indices |
| --- | ---: | --- | ---: | --- |
| High-sensitivity | 3.091261283 | 98 groups | 4 | 1 |
| Satellite v2 | 6.857877935 | 98 groups;1835 native pairs | 7 | 1,3 |
| OrthoFinder full | 9.947630559 | 99 checkpoint groups;1834 native pairs | 10 | none |

All whole-command screens pass, but that does not override interval flags.
All inputs contain645 proteins. These are small-fixture engineering results,
not accuracy estimates or scientific runtime comparisons.

## What The New Counters Establish

High-sensitivity interval1 has0.273399 unassigned CPU seconds under the
unchanged screen. Satellite intervals1 and3 have0.332489 and0.424304 CPU
seconds respectively. The corresponding average-core residuals exceed the
unchanged0.25 threshold. No run or flag was discarded.

For satellite interval3, host-minus-job is0.395465 CPU seconds, while the
separately sampled frontier measures only0.001811 outside-target CPU seconds
and retains a signed root-minus-frontier residual of0.377656 CPU seconds.
Across all three methods, outside-target frontier interval values range
from0.000799 to0.004611 CPU seconds. These observations do not identify
substantial CPU use in the sampled outside-job scopes as the explanation
for the flagged residual. They also do not prove the absence of interference:
counter windows differ, ancestor-direct tasks are unresolved, and transient
groups may appear and disappear between observations. Negative residuals
are retained too. No timing correction or causal attribution is justified.

Batch-step CPU per observation interval ranges0.028199-0.059998 seconds in
this array. Additional collector work is visible, but native wall-time
overhead has not been isolated by a matched observer-on/off experiment.

## Next Gate

Do not repeat this fixture selectively until it passes. The collector now
works end-to-end, but native observer overhead, non-CPU isolation, and a
defensible prospective inclusion policy remain unresolved. Existing
scientific scaling measurements remain unadmitted. Any policy revision
must be justified and frozen before new scientific measurements, preserving
the original screens and all earlier outcomes; small outside-job counters
alone cannot turn the current failures into accepted scientific timings.
