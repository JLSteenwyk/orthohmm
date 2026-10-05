# Native Factorial Reporting Result

The standard-library reporting CLI actually executes successfully on the
checksum-bound amended thirteen-identity plan, terminal reviews22427/22428,
failed-wrapper recovery22427 and separate frozen OrthoBench score22428.
See the [generated table](native_factorial_progress_20261004_v2/report.md),
[JSON](native_factorial_progress_20261004_v2/report.json),
[TSV](native_factorial_progress_20261004_v2/rows.tsv) and
[workflow](../NATIVE_FACTORIAL_REPORTING_20261004.md).

| Cell | Outcome | OrthoBench F1 (%) | Wall Seconds | Peak GiB |
| --- | --- | ---: | ---: | ---: |
| p0_c0_r0 | Failed wrapper; scientific outputs recovered | 69.7634 | 2289.7647 | 4.6140 |
| p0_c0_r1 | Native success; separately scored | 72.7050 | 3254.9227 | 4.6299 |

These are observed shared-host timings with unknown, potentially
tool-dependent distortion from competing analyses. The failed-wrapper row is
not a clean-success timing. CPU includes the wrapper bracket; lifetime peak
includes the native-step launcher, not pure algorithm RSS. Native intervals
exclude preparation, conversion and scoring. No timing correction, component
overhead, speed ranking, pooled median or independent accuracy is established.

An independent standard-library readback, without importing the exporter,
checks all13 JSON/TSV rows, both original terminal resource/score joins, foreign
CPU observations and all seven direct report/source pins. It confirms eleven
absent terminal reviews have null endpoints, not zeros or old cached scores.
The failed-wrapper recovery lacks an explicit reference-coverage numerator;
its coverage stays null. Successful22428 explicitly covers1,944/1,944 reference
genes. All920 helper sources bound to the live attempt remain unchanged.

The [joined regression run](native_factorial_progress_tests_20261004.xml)
passes465 tests, no failures/errors/skips,8.13 seconds. This includes41 new
reporter tests plus native adapter/execution/history/review/output/scoring/
recovery/QfO conversion/assessment tests. Tests and direct report readback do
not substitute for terminal native completion or raw scientific admission.

Fresh Slurm evidence while reporting confirms the same22429 attempt RUNNING
at20:52, with initial search52.78% complete. It is not scored and has no
supplied terminal review in this snapshot. That table state is not a claim
that the job is absent or terminal. Preserve the same job; review it and score
its output only after actual terminal evidence. No next identity is released.

The previously associated three configurations and original cached stage
costs remain separate. All broader publication requirements remain active;
this partial reporting snapshot is neither complete factorial evidence nor
a submission-ready manuscript/archive update. Existing mainPDF/rc4 and all
scientific settings, results, raw receipts and unrelated work are untouched.
