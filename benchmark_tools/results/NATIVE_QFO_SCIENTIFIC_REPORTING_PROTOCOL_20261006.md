# Native QfO Scientific Reporting

`benchmark_tools.export_native_qfo_scientific_scores` accepts separately supplied
ordinary native admissions and explicit recovered-science admissions. This is
an extension outside the frozen920-helper set, not a modification of the
success-only exporter or a synthetic successful measurement receipt.

For ordinary admissions it calls the existing strict reporter. Recovered
admissions must independently establish accuracy, match original conversion,
request/recovery/plan/native index, successful conversion/scoring allocations,
failed native scheduler outcome, execution/preflight/source/FAS bindings and
native endpoint arithmetic. Resources stay null; timing admitted/eligible false.
Report rows identify accuracy status and measurement status separately. Missing
cells remain null/blank/Unavailable, not zeros or inferred running states.
Reject duplicated cells even across normal/recovered lists, cached participants,
reused P1C0R0, source/reference changes and failure/success relabeling.

The table retains native precision and recall, assessed relation counts,
endpoint-specific error fields, all-input relation coverage, resolved versus
clique prediction semantics, FAS protocol/sample limits and original provenance.
Only VGNC/SwissTrees/TreeFam-A summaries are F1; GO/EC similarity and FAS are
different statistics. The six-metric mean is project-defined and secondary,
not official QfO F1. No native SEM becomes a paired-family confidence interval.

This performs direct reporting/source/arithmetic checks, not another transitive
raw admission, native inference, score execution, bootstrap or biological
validation. A successful terminal-review label is not a substitute resource
comparison. Shared-host effects remain unknown/tool-dependent; no timing
correction or isolated-efficiency ranking. Output must be a fresh directory.

316 joined tests pass3.48s, including60 new reporter cases and existing native
reporting/recovered conversion/assessment/recovery/endpoint contracts. Covers
mixed normal/recovered admissions, timing-null flags, precision/recall/non-F1
semantics, incomplete/changed/nonfinite arithmetic, input denominator, pair
coverage bounds, duplicate cells, receipt hashes, missing values and no overwrite.
The separate reproduction-guide contract suite passes12 tests in1.37s.
No real recovered-score table is produced at this preparation checkpoint:
actual scoring22448 remains live and independent validator22449 waits afterany.
Use only actual independently admitted results after terminal completion.

The [actual scorer workflow](RECOVERED_NATIVE_QFO_ASSESSMENT_PROTOCOL_20261006.md)
was committed/pushed9bdeb160 before launch. Its fixed manifest retains the
Nextflow22.10.8/workflow19.10.0 compatibility warning. Initial live-start observer
looks for the shell token in a Java-child argv and fails; corrected observation
checks `nextflow.cli.Launcher` plus the exact frozen command suffix. The failed
assertion is recorded, not an analysis restart or a compatibility certification.
