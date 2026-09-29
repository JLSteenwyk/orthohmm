# Deferred Independent Completion Audits

## Prospective Scope

Submit one bookkeeping/audit job after both 22377 and 22378 become terminal,
using `afterany:22377:22378`. This does not repeat either native workload,
retry an inference failure, authorize the timing panel or change any scientific
configuration. Both retained audit commands are attempted independently once;
one failed audit does not suppress the other. Failure, timeout and mismatch
outputs remain retained, with no automatic replacement attempt.

The [machine-readable protocol](publication_completion_audits_protocol_20260929.json)
is frozen before submission at SHA-256
`8dd65820a2d9ca1a40dd8166391f0c8f3240358e0e76ac951657bc775c7c589e`.
It binds the existing 822 unchanged calibration harness files plus the new
orchestrator, its exact trusted audit interpreter, submission script and
commands. It does not upgrade the native/scientific runtime or grant general
security clearance to the audit environment. Sources are rechecked before
and after the audits; do not silently refresh pins during this attempt.

Resources: `bizon`, four CPUs per task, 32 GiB RAM, six-hour allocation ceiling,
no requeue. Both audit subprocesses have a 7,200-second timeout. Single-thread
numeric-library settings prevent oversubscribing this audit allocation. These
resources are not inference measurements and are not pooled into timing tables.

## Stages And Acceptance

1. Run the unchanged `admit_restored_archive_ob` against job 22377's retained
   directory and fresh `independent_admission` output. It verifies execution
   and installed packages, invokes independent readers, recomputes all RefOG
   score objects, and compares full partitions and four native TSVs against
   admitted 22376. Audit completion and reproduction equality are separate.
2. Run unchanged `audit_threadripper_observer` against job 22378's retained
   directory, the externally pinned calibration protocol, and fresh
   `independent_audit.json`. It replays raw evidence and evaluates the previously
   frozen accounting/cadence checks, without admitting scientific timing.

Capture parent allocation/batch/step accounting first. For each stage, success
requires zero CLI exit, a matching parent job ID and recognized result status,
the appropriate semantic check, and completed zero-exit parent accounting.
OrthoBench additionally requires `reproduction_equal: true`; CLI exit 0 alone
is insufficient. Calibration additionally requires `calibration_checks_passed`.
Failed/incomplete parent accounting cannot produce success. Both outputs and
logs are retained regardless of the combined outcome.

The runner refuses an existing control directory or audit result. Existing
partial audit directories are left intact and cause the original validators
to refuse repetition. Do not invoke manual versions of these commands while
the deferred audit is pending/running. After it finishes, inspect its results
rather than executing the same audit again.

Owned timeout/interruption cleanup signals only the launched audit process
group and reaps its subprocess. It does not signal either original job or
unrelated services. No scheduler configuration or DGX operation is performed.

## Retained Outputs

Control directory:
`benchmarks/work/publication_completion_audits_20260929`

- `parent_accounting.json`: exact query and observed scheduler output.
- `restored_orthobench.log` / `.json` and `observer_calibration.log` / `.json`:
  commands, exit/failure/timeout, semantic and parent-accounting checks, and
  result/log file identities.
- `complete.json`: protocol/source identities, scheduler job ID and both
  stage outcomes. `completion_audits_passed` is a bookkeeping result, not
  publication readiness or controlled timing admission.

Scientific output remains under the original run's `independent_admission`;
calibration output remains under its original run's `independent_audit.json`
and independent scoped-resource receipt. Large native inputs/results are not
copied into Git or into this audit job's accounting summary.

## Validation Before Submission

74 focused orchestrator, existing scientific-admission and calibration tests
pass. Coverage includes retained mismatch despite exit 0, wrong result jobs,
nonzero exits, failed calibration, parent allocation/batch/step failure or
missing/malformed accounting, fresh-output refusal, owned-child timeout/reaping,
source/protocol binding and both stages attempted despite a failed first stage.
Orchestrator composition uses explicitly synthetic audit commands/accounting;
it is not evidence that the two pending native jobs completed or passed.
Shell syntax validation passes. The full publication goal remains active.
