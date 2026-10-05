# OrthoBench Terminal Followup Queued As22432

After pushing17d24eea's tested prospective
[worker/protocol](../NATIVE_ORTHOBENCH_FOLLOWUP_20261005.md), actually submit
the review/score job on hold. Its parent is the same live native22431,
index4/P1C0R1, not a new inference or retry. No source/runtime/scientific
parameters of the parent are changed.

The [held submission](native_orthobench_followup_submission_22432.json)
records the actual argument vector, source commit, returned job ID, worker,
batch and original request pins. An empty owned followup-name queue is
checked before submission, preventing an uninspected duplicate. The
[stored Slurm batch script](native_orthobench_followup_batch_22432.sh) is
byte-identical to the tested765-byte source script. Worker SHA256 is
`7a5962eea9d775b9143e700e8afaf7952c6f88ff76e9f96030ff900c6a790303`;
the original request remains
`436541a8977289cf190c1d3918eea1930e5eb7c5ef518421670971ae7dfc583f`.

The [independent pre-release check](native_orthobench_followup_pre_release_22432.json)
freshly verifies owned PENDING/JobHeldUser, afterany22431(unfulfilled),
the bizon/gpu request, twoCPUs/8GiB/two hours, no requeue/restarts/arrays,
command/cwd/comment, absent review/score destinations, source/request and
generated-script bytes, every submitted argument and all920 unchanged
native helpers. AvailableRAM662,367,207,424bytes. This rechecks submitted
argv from the recorded submission; the controller reports the script path,
not an independently retrieved full executed argv. No full runtime closure
claim follows.

An optional JSON-controller client probe fails with missing serializer
plugin/exit139. Fresh text-controller observations work and are used instead.
The native and followup jobs are unaffected; no job is restarted and no
Slurm configuration/plugin/service is modified to repair this optional probe.

Only after those checks, actually release22432 once. The
[post-release observation](native_orthobench_followup_released_22432.json)
confirms pending, no longer user-held, with the same unfulfilled dependency,
owner and worker comment. Parent22431 remains RUNNING with livePID2805435
and affinity0..31. The pre-release receipt deliberately remains historical
and says release not executed; it is not overwritten with later state.

The followup's original pinned local outputs are:

- `benchmarks/work/native_factorial_terminal_review_22431/review.json`
- `benchmarks/work/native_factorial_terminal_review_22431/followup.json`
- `benchmarks/work/native_factorial_orthobench_score_22431/score.json`

These are expected destinations, not current successful results. On resumption,
inspect both actual handles22431/22432. Never manually duplicate the review
while the dependent worker is pending/running. Wait for actual whole-job
terminal evidence, inspect the full native review/separate score and worker
outcome, then publish/read back results before preparing index5's complete
reviewed-prefix and held-release checks. If a review/scorer fails, retain it
and diagnose that failure without restarting native inference or silently
retrying into an existing destination.

94 joined tests and matching-runtime import/bash syntax are preparation
evidence, not a completed real review. No fresh native score or next-identity
release is admitted. Postprocessing resources are separate from native
inference; shared-host contention has unknown potentially tool-dependent
timing effects. The full publication goal and remaining scientific/native/
distribution requirements stay active and unproven.
