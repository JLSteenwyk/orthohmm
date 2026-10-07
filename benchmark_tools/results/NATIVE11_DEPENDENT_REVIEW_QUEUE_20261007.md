# Native11 Dependent Review Queue

The active execution contract already authorizes submitting, inspecting,
dependency-chaining and releasing this project's own jobs. This prospective
workflow change replaces the previous plan to defer reviewer submission until
observing terminal native11 state. It does NOT permit reviewer execution on
live/partial native outputs, scoring failed outputs or automatic retries.

Queue the unchanged, bash-checked
`allocated_native_review_23985_20261007_v1.sh`, SHA256
`8093add87b49fb4aa470728d084e925b00d595d9fff09395cb0ff0ffb32462a0`,
with `--dependency=afterany:23985`, initially held. The native request remains
`7bf63b80bd5932b9edbd1b2c5ff3fb77f5557e4f6c64e077045d6a50c8d366a1`.
Both script and original reviewer are unchanged. No live scientific source,
native job allocation, unrelated job or scheduler configuration is modified.

Before submission, inspect the existing native job, own reviewer job-name
handles, fresh output namespace and script/request hashes. Submit once with
two CPUs,128GiB,six hours,no requeue, distinct stdout/stderr names and the
source digest in the scheduler Comment. Inspect the held job's owner, command,
directory, allocation, dependency and requeue policy before one release.
Record actual observations/submission/release in separate receipts; preserve
any failed or ambiguous outcome without blind resubmission. The script has
its existing fresh-at-execution memory-capacity and distinct-job checks.

`afterany` permits classification of unsuccessful native terminal outcomes
as well as successful ones. It does not mark either outcome scientifically
admissible. The unchanged reviewer requires the actual native terminal state,
validates the complete retained history/request/environment/accounting and
handles failed-output evidence separately. Existing reviews/namespaces must
not be overwritten. The scheduler dependency prevents early execution; the
reviewer's own terminal validation remains an independent execution gate.

After releasing, inspect the same new handle and native23985. Pending due to
Dependency is a verified waiting state, not a blocker or a request for another
user resume. At actual reviewer completion inspect its report/failure and
original next-identity authorization. Conversion/scoring/admission or genuinely
unrun native12 remain gated by their existing contracts. No automatic downstream
launch or claim of success follows merely from queueing this review.
