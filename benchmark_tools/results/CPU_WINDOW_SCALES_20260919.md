# CPU Accounting Across Window Scales

Recomputed counters directly between retained endpoints at fixed5-,10- and
30-sample widths, starting from point0. Kept every final partial window and
all original one-sample flags. Did not sum overlapping residuals or choose
window alignment using flagged locations. This is post-outcome engineering
description, not a new timing-admission rule.

| Native method | One-sample narrow flags | Five-sample | Ten-sample | Thirty-sample | Whole-observation residual cores |
| --- | ---: | ---: | ---: | ---: | ---: |
| OrthoHMM high sensitivity | 0/551 | 0/111 | 0/56 | 0/19 | 0.068205 |
| OrthoHMM satellite_v2 | 21/819 | 1/164 | 0/82 | 0/28 | 0.100014 |
| OrthoFinder full | 1/622 | 0/125 | 0/63 | 0/21 | 0.059319 |

Satellite_v2 maximum residuals on the three longer grids are0.270538,
0.242144 and0.224560cores. All native whole-observation narrow screens
pass their unchanged numerical checks. Whole-observation windows bracket
wrapper work; they are not exact native command boundaries. Exact-command
original screens are separately retained in the report.

All six native-only controls remain unflagged at every tested scale. Each
known-contended control flags all5 five-sample windows, all3 ten-sample
windows and its single thirty-sample partial window. Control commands contain
only21 one-sample intervals, so their thirty-sample result is a whole-run
partial window, not a30-second exposure. Their whole-observation residuals
are0.950954,0.958079 and0.957512cores.

Longer aggregation suppresses the historical native flags while retaining
these sustained positive controls. It can suppress brief real interference
as well as non-atomic counter noise, so these results cannot identify the
historical cause or retroactively accept timings. Unchanged screen formulas
at longer durations are diagnostic only, not newly calibrated gates.

Raw report: `benchmarks/work/cpu_window_scales_20260919.json`, SHA-256
`34eb1c0e88149664a90a7ee886577a73e46e5c22d941384d0262f11c0ddd2b0d`.
Retained gzip: `cpu_window_scales_20260919.json.gz`, SHA-256
`c8beb10d1d9f749bafc030a4bab505c37ec2aee24854fb2b76faefab6ea3391b`.
Input audits are checksum-pinned; original native measurement hashes are
checked before/after reading. Seventeen focused scale/bracket tests pass.

Separately froze `DUAL_COLLECTOR_OVERHEAD_PROTOCOL_20260919.md` and derived
the complete fresh18-task overhead plan without changing native workloads
or the prior5% median/10% pair budgets. Protocol SHA-256:
`7df6f18476d20d9e96ab1a2f1fd69956e42347bf519ba93f35c0edeb1c4e761a`.
Plan SHA-256:
`1160a669a4033c66e7bdde5baddf429e83b9328d3821f99291fd29210cac8ab9`.
An additional test verifies actual parent-plan commands/order/resources are
preserved after path relocation. No new overhead tasks have been submitted;
pinned launcher/deployment remains next. No scientific timing admitted.
