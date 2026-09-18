# Corrected OrthoMCL QfO Assessment

## Protocol

OrthoMCL is added to the existing corrected comparator workflow without
changing the six frozen QfO endpoints, reference data, scoring containers,
command construction or endpoint interpretation. Its participant is
`qfo_corrected_orthomcl`; its dedicated native work directory is
`qfo_benchmark/w/qc_mc`, and output is
`qfo_benchmark/scoring/corrected_orthomcl`.

The scoring executor is frozen before submission at
`benchmarks/work/publication_qfo_corrected_orthomcl_assessment_v1`; its final
revision is recorded with the submission below.
Converter: `b213c3c1087510e7e965afa14ffca35d14b8ca1a` at
`benchmarks/work/publication_qfo_corrected_orthomcl_pairs_v1`.
The already queued OrthoFinder, FastOMA and other comparator executors are
not modified, and their existing provenance formats are preserved.

Scoring requires completed 2-CPU/64-GiB conversion on bizon, matching job and
source identity, the frozen converter revision, positive pair count with
zero mapping loss, byte-identical raw and filtered pair files, matching
group and ungrouped counts, and final-group clique semantics. An empty pair
file is explicitly rejected by this assessment workflow, not silently
treated as a successful scored method.

The conversion's retained group audit is rechecked against the source native
admission. Its content, paths, hashes and query-failure diagnostics must
match; its checked source and input records are incorporated into scoring
provenance. The independent score admission reconstructs the same records
and commands, validates terminal scoring and conversion accounting, checks
the immutable output inventory and native task trace, and validates all
endpoint output files with the existing QfO validator.

The admitted summary explicitly retains OrthoMCL pair semantics, group
coverage, source query diagnostics and the group-audit record. This is
score/provenance validation, not independent biological validation or proof
that absent BLAST hits were repaired. The six-endpoint mean remains a
project-defined secondary summary, not an official QfO F1 score.

The corrected-table exporter now has an explicit `orthomcl_1_4` adapter
which preserves these semantics and coverage/diagnostic fields in its
machine-readable row. It rejects inconsistent bindings or mapping loss.
Pending methods remain missing rows, never zero scores.

## Verification

94 focused assessment/admission tests pass, including all six comparator
identities, OrthoMCL corruption cases, additional native-audit evidence,
pending-job rejection and mocked complete score-admission orchestration.
Randomized/native-fixture pair correctness is covered separately by
QFO_CORRECTED_ORTHOMCL_PAIRS_20260918.md. Mocked metric validation is not a
production result.

The first full-suite run exposed a missing OrthoMCL corrected-table adapter
(one failure,4,033 passes). The adapter was implemented with six additional
tests rather than weakening the pipeline/exporter completeness check.

An actual check-only scorer invocation against pending conversion21753
rejected missing terminal accounting before reading partial conversion
output. Both new batch scripts pass `bash -n`. Full-suite verification and
submission identities are recorded below once available.

## Scheduling

- Assessment batch: `qfo_corrected_orthomcl_assessment_batch_20260918.sh`,
  8 CPUs/64 GiB/24 hours on bizon, no requeue, `afterok:21753`.
- Independent admission batch:
  `qfo_corrected_orthomcl_score_admit_batch_20260918.sh`, 2 CPUs/64 GiB/
  24 hours on bizon, no requeue, `afterany` of the assessment job. The larger
  verification allowance accommodates rehashing retained large native
  provenance; it does not change the scoring endpoint or resource settings.

The pair-manifest checksum is computed after dependency satisfaction.
Expected final admission is
`benchmarks/work/qfo_corrected_orthomcl_assessment_admission_20260918.json`.
No corrected OrthoMCL score or successful production assessment is claimed
by preparation or queue submission.

## Submitted Jobs

Full suite after adding the exporter: **4,040 passed in87.10s**, with legacy
native runtime tests enabled. After pinning the final scoring revision,
120 focused assessment/admission/export tests pass. An actual invocation of
the frozen admission executor correctly rejects pending21754 before reading
or creating score artifacts.

- **21754**, submitted scheduler time `2026-09-18T10:01:41`,
  `afterok:21753`,8CPUs/64GiB/24h/bizon, no requeue.
  Frozen scoring revision `5efb206b23a44f85386d8cb7e90f3e815a3d162a`.
- **21755**, submitted scheduler time `2026-09-18T10:01:52`,
  `afterany:21754`,2CPUs/64GiB/24h/bizon, no requeue.
  Frozen validator revision `4af373076beb138a086c159c7950e411762197d1`, at
  `benchmarks/work/publication_qfo_corrected_orthomcl_score_admission_v1`.

Both dependencies and resources were confirmed by `scontrol`. The unused
provisional scorer checkout at e1b4944 was advanced, while clean and before
submission, to the fully tested5efb206 snapshot including the table adapter.
No running executor or previous method's frozen worktree was changed.
