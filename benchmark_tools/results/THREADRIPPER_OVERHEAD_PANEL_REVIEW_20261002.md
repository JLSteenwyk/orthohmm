# Conditional Engineering Panel Review

Added `benchmark_tools.review_threadripper_overhead` to combine independent
native/resource/output replay with the four existing recorded review categories.
It does not submit jobs or certify the conclusions of runtime/environment
reviews. No real overhead-budget result is available.

## Contract

Inputs are externally pinned plan, attempt history and review history, plus
the existing baseline/launcher/scheduler bindings when attempts exist.
The plan must preserve the fixed private scientific parent. Raw auditing is
recomputed with the existing auditor; its unrun/failed outcomes, output
identities, signed ratios and complete-cell arithmetic remain intact.

Review history uses `threadripper_overhead_reviews_v1` with exact `plan`
and `attempt_history` records and one ordered `sessions` reference for
every recorded attempt, not only successful or favorable ones. An empty
session list is valid only for an empty attempt history.
Each engineering session/review retains the existing pair/arm/job/plan schemas.
The controller stdout must equal that attempt's pinned raw scheduler record.
All four review categories require explicit nonblank `review_reference`.
The output category must pin its native audit, whose outcome, checked outputs
and wall time must agree with the independently recomputed replay.

A conditional engineering decision requires all **54 successful, reviewed,
output-consistent tasks**, a completely reviewed contiguous prefix and the
existing numerical budget result. It preserves the <=0.10 per-pair and
<=0.05 complete three-repeat median limits. Missing or ineligible evidence
produces `engineering_budget_passed: null`, not zero, partial medians or
a passing/failing complete-panel budget. A complete reviewed panel can report
either a numerical pass or a failure; slow results are not discarded.

The `reviewed_runtime_environment_complete` flag means all recorded review
categories are bound and consistent with replay, **not** independently certified
environment/runtime validity. `review_conclusions_independently_certified`,
scientific timing admission, publication readiness, automatic retry and
next submission stay false. The raw auditor's separate admission flags remain
unchanged. Authenticity and scientific validity of external review conclusions
remain review responsibilities; this helper is not a quiet-host oracle.

Invoke `python -B -m benchmark_tools.review_threadripper_overhead` with
`--plan`/`--plan-sha256`, `--attempts`/`--attempts-sha256`,
`--reviews`/`--reviews-sha256` and a fresh `--output`.
For attempted work also supply `--baseline`/`--baseline-sha256`,
`--launcher`, `--scheduler-command` and `--scheduler-cwd`.
This is an audit command, not a launch command.

## Executed Evidence

Source is committed at `e617c4061b1d8d833b8292d1cac84b5b17b5163d`.
Four-module regression passes **137 cases in 24.16 seconds**, zero failures,
errors or skips. The complete pass/fail cases use fabricated review records
and a test-double raw auditor; they are **not actual 54-task measurements**.
Real unrun/CLI tests use copied, hash-checked local evidence, not workstation
paths. Mixed plans/history, scheduler discrepancies, native audit differences,
unsupported references, production schemas and missing/failed reviews cannot
produce a passing decision.

Run the actual CLI under isolated test Python on the unchanged retained private
plan/empty attempt history and a new
[empty review inventory](threadripper_native_overhead_reviews_unrun_20261002.json).
It exits zero with `engineering_panel_incomplete_or_ineligible`: all 54
tasks remain unrun, 27 pairs and nine cells have no complete medians, and the
budget is null. No actual allocation, host/environment observation or inference
executes. Zero CLI exit is successful audit execution, not engineering approval.

Separate stdlib readback verifies 14 unique evidence files, raw/combined task
counts and missing decisions, flags, process/JUnit data and source bytes against
Git. Preserve a one-line readback SyntaxError before execution; no output file
or reviewer rerun occurred. See the
[compact validation receipt](threadripper_overhead_panel_review_validation_20261002.json):
3,897 bytes, SHA256
`7ffb9e91f95e8f04da3ad082840fde5d730ae09285bbf73d05ed7acdc75248b9`.

## Outstanding

This closes review combination, not actual overhead or production timing.
Final source/runtime binding must include the new helper; old startup inventories
are not current approval. Real environmental handoff/policy/readiness, isolation,
all 54 engineering outcomes and the separate 27 production runs remain required.
No DGX, unrelated-job/service change, shared upgrade, scientific retuning,
archive/PDF rebuild or timing launch occurs.
