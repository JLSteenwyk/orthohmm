# Candidate Import-Order Correction

## Preserved Failure

New job 22384 used committed executor `3ac5e52395d51df93246fd312e5c112e5efbfdfd`
and the explicitly selected admitted seed. It is terminal FAILED, exit 1:0,
2 CPUs/64 GiB on bizon, scheduler interval 20 seconds. The traceback records
`Wrong frozen candidate scientific import`; no output root was created and
the pinned source raises before the candidate builder can run. The failed
attempt is not a candidate or accuracy result.

The [submission receipt](qfo_cpm_helper_candidates_submission_22384.json), SHA256
`bc9e5f919b2c68ad75a6246a365d8b19868a1c589462e3662a352d8b6a438d93`,
and [failure receipt](qfo_cpm_helper_candidates_failure_22384.json), SHA256
`4ef263a9b34308c2b4c0f73ff410f650342e4348b89a74bd52bf634ae239ea70`,
retain executor/source/script/input identities, scheduler state, log, GNU time
and a new bounded fresh-process import-order diagnostic. Its diagnostic is a
reproduction, not an original worker module snapshot. It shows the status helper
`run_blast_recovery_batch` transitively caching the executor's OrthoHMM package
before the frozen launcher is selected; later engine/accuracy imports therefore
resolve from the executor, and the guard correctly rejects them. This specific
engineering failure is unrelated to the historical native-refinement crashes.

## Prospective Corrected Attempt

Move only the status-helper import after successful frozen engine/accuracy
selection and before checkpoint auditing. Retain the unchanged status writer,
frozen scientific-origin guard, numeric auditor, builder, settings and runtime
comparisons. Do not delete/reload cached scientific modules or accept a wrong
package. Two fresh-process real-import regressions exercise actual transitive
helper imports: a clean selected package succeeds, a previously imported wrong
package remains rejected. Their scientific execution and numeric data are
fixtures, not biological accuracy or equivalence evidence.

This deterministic source correction justifies one separately identified new
candidate attempt, not an automatic retry of unchanged code. First commit/push
the fix, tests, failure receipts and this amendment. Create a new detached
executor at that exact revision; use the same hash-bound prospective candidate
protocol, readback and committed launch script. Record this correction document
as well as the new executor/script/source identities in its submission receipt.
Original output root is still absent; require it to remain absent at launch.

Keep original Python, frozen QfO core/launcher/native libraries, numeric hits,
CPM 0.12/seed 4 and fixed satellite_v2 candidate parameters. Keep 2 CPUs, 64 GiB,
4-hour limit, local bizon, no GPU/requeue, all one-thread library settings.
No reference labels, new search, optimizer/refinement, parameter tuning or
timing-panel change. Preserve failure 22384 and all older failed/cancelled jobs.
Stop and retain any new failure; do not automatically retry. Success remains
pending independent candidate admission and all downstream accuracy gates.

GNU time's 19.17-second failed preflight interval is descriptive shared-host
history, not controlled comparative or end-to-end runtime. Timing remains
deferred, as do the other publication requirements.
