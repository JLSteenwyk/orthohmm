# Terminal Native Ablation Resources Consolidated

The [new table](native_factorial_terminal_resources_20261009_v1/resources.md),
[JSON](native_factorial_terminal_resources_20261009_v1/report.json) and
[TSV](native_factorial_terminal_resources_20261009_v1/resources.tsv) join all
13 retained native attempt identities. Eleven have resource observations;
ten have admitted accuracy. These are distinct counts, not complete successes.
The table does not manufacture scores, replace failed attempts or measure new
timings. The original v7 resource snapshot and every review remain unchanged.

| Later Identity | Native Outcome | Resources | Accuracy |
| --- | --- | --- | --- |
| 7, P0/C0/R1 | Measurement failure; scientific outputs recovered | Unavailable | Admitted |
| 8, P0/C1/R0 | Native success | Reviewed native-command observation | Admitted |
| 9, P0/C1/R1 | Pre-native CPU-binding failure | Unavailable | Missing |
| 10, P1/C0/R1 | Allocated native success | Reviewed native-command observation | Admitted |
| 11, P1/C1/R0 | Native success; subsequent scoring OOM | Reviewed native-command observation | Missing |
| 12, P1/C1/R1 | Native SIGSEGV | Failed-command observation only | Missing |

The six fresh OrthoBench rows and QfO index 6 keep every old resource number.
OrthoBench index 0 retains its failed-wrapper observation and recovered accuracy;
it is not relabeled clean-success timing. Later review schemas remain explicit
and distinct. The reporter does not turn their false scientific-timing-admission
flags into a successful timing admission.

## Actual Execution And Readback

Prospective source/tests committed/pushed at
`9598abc900556ac94d750e557e9e587ad9f63367` before actual export. The final focused
and adjacent test suite passed 111 cases in 2.11s. First invocation failed 13
new cases because the measurement-failure schema omits `resource_scopes`, rather
than storing an explicit null. Corrected only the new missing-resource consumer
and retained the original failed measurement. No production output preceded
the tested correction.

Production export executed once with isolated `-I -S -B` Python 3.10 and
sanitized environment, returning zero. `report.json` is 16135 bytes, SHA256
`061f299cbeb2a1f2eeb69410ff548f4a90d5e9537767545cf02b23ba2a86af65`.
The [actual execution/readback receipt](native_factorial_terminal_resource_execution_20261009_v1.json)
records the command, output and independent standard-library readback code.
Without importing the exporter, that readback rehashes all eight direct inputs,
checks all 13 resource vectors against original reports/reviews, all accuracy
flags against the current QfO snapshot, all 156 TSV cells, every displayed
number and unchanged baseline rows. It exits zero. This is independent
projection/arithmetic verification, not a fresh raw accounting replay.

## Scientific Scope

Native intervals exclude preparation, conversion and scoring. CPU is the task-
subtree bracket including its wrapper; peak is native-step lifetime memory
including the launcher. Failed-command observations do not estimate successful
inference costs. No pool, corrected timing, component-overhead claim or isolated
speed rank is produced. Shared-host distortion remains unknown and potentially
method dependent; no DGX or quiet-window gate is reinstated.

The [six associated historical OrthoBench observations](factorial_native_resource_linkage_20261004/linkage.json)
and [historical QfO high-sensitivity cost](qfo_native_configuration_cost_20261004/association.json)
remain separate. They cover the other three configuration identities, with
different repeat/runtime/memory scopes and partition checks. None supplies the
cost of an original cached factorial execution; all 16 unavailable original
cached full costs stay unavailable. The current table closes the stale terminal
join identified in requirement 3.4, not every broader cost/provenance limitation.

The rc6 archive and current frozen manuscript are not rebuilt to absorb this
addendum. No method/default or primary endpoint changed. Appropriate paired
uncertainty and biological-stratum limitations remain open; publication readiness
is false. Next inspect whether original comparison protocols and exact common
SwissTrees family records support native-cell versus OrthoFinder uncertainty,
without selecting a favorable cell or silently expanding corrected endpoints.
