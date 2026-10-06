# Native Reviewer Crash: Diagnostic Launched

## Original Outcomes

Fresh accounting confirms native22444 COMPLETED0:0 in11:42:59; its native
step completed0:0 in11:18:09. Original reviewer22445 FAILED0:11 in31:25.
Original downstream22450/22451/22452 FAILED1:0 at their predecessor-accounting
gates, before conversion, scoring or admission. Preserve their logs and gate
results; do not claim the native result is scored or terminal-reviewed.

[Actual failure observation](native_review_sigsegv_observation_20261006.json)
binds the original request/plan, all five partial reviewer files, eight original
stdout/stderr logs and three refused downstream gate results. The kernel
records Python PID818632 segfault15:20:23; retained-binary addr2line identifies
_Py_INCREF, not the cause. No final reviewer success/failure JSON exists.
Partial runtime/resource/environment reports are not a complete admission.
Observation SHA256:
bc8c743afaea219c9c2b1a94a3f6f3fb44919e56fdcac16c1680550d453bb481.

## Actual Diagnostic Launch

[Protocol](NATIVE_REVIEW_SIGSEGV_DIAGNOSTIC_PROTOCOL_20261006.md) and
[batch](native_review_sigsegv_diagnostic_20261006.sh) committed/pushed5200f437
before submission. One new job22734,2CPUs/32GiB/6hours/no-requeue/bizon/gpu,
runs the unchanged original semantic validation CLI in the retained Python3.10
review venv with fault reporting enabled. No native inference/resource replay,
conversion/scoring/admission is rerun by this diagnostic.

Initial held-inspection inline observer refuses Slurm's NumNodes=1-1
representation because it expects literal1. It releases nothing and writes
no receipt. Re-poll and inspect the SAME held22734, accepting only equivalent
one-node representations with exact node=1 TRES; do not resubmit. Actual
[held receipt](native_review_sigsegv_submission_22734.json) records this
observer failure and the fresh validated identity/resource envelope. Submission
SHA256:4f1df72adca42257ec48d8b3585a4e1f4f462d6aefb46098fd913df5c9db4d14.
Initial submission stdout is not retained as a file; do not reconstruct it.

Observation/submission milestone0804577d pushed before one release.
[Release receipt](native_review_sigsegv_release_22734.json) records before,
release command and immediate PENDING/Reason=None observation. Fresh accounting
subsequently confirms RUNNING; owned PythonPID1198962 is R, with56seconds CPU
at56seconds elapsed and about999MiB RSS. This is a verified live diagnostic,
not a successful outcome. Do not read unfinished result JSON as admission.
Release SHA256:f10cc8db4356ea5353eaa4343eab7894005fa1800b4ff23bc6d1cd53435bd28e.

All920 frozen helper sources, plan/evidence/request and original interpreter
bindings checked before release; no frozen files modified. Scientific3.10
imports Bio1.87/NumPy2.2.6/psutil7.2.2. Available RAM at failure observation
645,655,040,000bytes; current safe capacity rechecked before release. Nearly
full swap is recorded, not a reason to stop unrelated work or invent isolation.

## Validation And Next Action

Existing output, terminal-review and conversion contract tests:156passed3.04s,
zero errors/failures/skips, XML time2.986s. Retained
[test receipt](native_review_signal11_existing_contracts_20261006.xml).
Tests run in the retained Python3.12 test environment; scientific Python3.10
has no pytest and is not modified. Bash syntax passes. These tests do not
reproduce or explain the production crash or validate an unfinished diagnostic.

Poll original22734 until terminal. Inspect its actual fault stack/error or
validated output only after completion. Success would prove standalone
semantic validation in that invocation only; a separately documented full
postprocessing repair and all unchanged gates remain required. Keep native
index9 gated; no automatic retry, substitution of diagnostic admission or
rewriting of the original22445/22450/22451/22452 failures.

Elapsed times are shared-host diagnostic observations. Competing CPU,
memory-bandwidth and I/O demands may distort them by an unknown, potentially
tool-dependent amount; no corrected time, isolated ranking or inference cost.
The full publication goal is active, including remaining native cells,
interactions/matched search, appropriate uncertainty, independent validation,
provenance, TreeFam limitations and executable manuscript/archive/release.
Completion remains unproven.
