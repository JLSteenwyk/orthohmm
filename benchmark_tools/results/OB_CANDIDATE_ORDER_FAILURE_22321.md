# Candidate Diagnostic Startup Failure

Job 22321, submitted from commit `d0efe55`, terminated `FAILED`, exit `1:0`,
after one second. The [receipt](ob_candidate_order_failure_22321.json) retains
scheduler evidence, the log, original plan and a saved copy of the submitted
driver whose hash matches that plan. The
[submission](ob_candidate_order_score_submission_22321.json) is retained.

The driver imported the strict candidate auditor into the private scientific
runtime. That auditor imports a historical FASTA/scoring module requiring
Biopython, which is absent from this intentionally minimal runtime. Import
failed before `started.json`, hit loading, candidate expansion or any arm
output. This is a diagnostic harness failure, not an OrthoHMM inference
failure or biological result. No retry or resume was performed.

The correction keeps native execution separate from external validation.
The isolated driver writes candidate partitions, summaries and output hashes;
the existing strict seed/merge auditor must run afterward in the audit
environment. Biopython was not installed into or otherwise used to change the
scientific environment. Neither candidate algorithm nor scientific defaults
changed.

A generated four-gene fixture now executes all five arms through the same
isolated driver and passes all five external candidate readbacks. The fixture
includes reordered nonself hits, small score perturbations and one self hit.
The receipt pins its plan, execution, report and independent readback.
Fifty-six focused tests also pass. Fixture success is not full-data admission.

A revised full-data plan has been prepared but **not submitted**:
`benchmarks/work/ob_candidate_order_scores_v2_20260926/plan.json`, SHA-256
`cc3bdcef846c68b5d6a2e524a976d68325cfd17e184fdb2934eac618011a9ae8`.
The five scientific arms, resource limits and one-attempt rule remain as
specified in the [original protocol](OB_CANDIDATE_ORDER_SCORE_PROTOCOL_20260926.md).
Any new submission must use this distinct plan and preserve job 22321's failure;
the old output directory must not be resumed. Terminal scheduler success,
independent candidate checks and baseline comparisons remain required before
interpreting full-data results. No score/order mechanism has yet been tested
on the full OrthoBench input by this factorial.
