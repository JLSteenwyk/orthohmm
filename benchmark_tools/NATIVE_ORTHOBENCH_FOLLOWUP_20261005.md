# Dependent Native OrthoBench Terminal Followup

The same index4/P1C0R1 native job22431 is running. Prepare its already-required
terminal review and frozen scoring as a separate dependent postprocessing
job, without waiting for another interactive turn after inference finishes.
This changes execution scheduling only, not the 13 native identities,
scientific settings, endpoints, resource scopes or frozen helper sources.

Use [the small worker](finish_native_orthobench_attempt.py) and
[its batch envelope](finish_native_orthobench_attempt.sh). Submit on hold
with `afterany:22431`, on bizon/gpu, 2 CPUs, 8GiB, two hours, no requeue.
Before release, independently check the owned held job's exact dependency,
argv, source/request checksums, command, working directory and resource
envelope. Retain its generated Slurm batch script and submission/controller
observations. No scheduler configuration or unrelated job is changed.

The worker requires matching Python3.10 and explicit source/request hashes.
Both destinations must be fresh, direct, distinct repository paths. The
unchanged terminal reviewer freshly checks the whole parent job, request,
runtime, resource replay, shared-host observations and native outputs. Only
its bound native success reaches the unchanged separate scorer. Failed
native attempts remain unscored; unknown outcomes are refused, not relabelled.
Partial reviews and scorer errors remain retained. No automatic retry,
inference launch or next-identity release occurs.

For the first planned invocation, use original request
`benchmarks/work/native_factorial_launch_20261004/request_04_receipt_amended.json`,
SHA256 `436541a8977289cf190c1d3918eea1930e5eb7c5ef518421670971ae7dfc583f`.
Fresh destinations are
`benchmarks/work/native_factorial_terminal_review_22431` and
`benchmarks/work/native_factorial_orthobench_score_22431`.
The actual submission receipt must pin the final pushed worker source hash;
do not substitute a future edited worker.

`afterany` permits inspection of failed parent jobs; it does not make them
successful. Pending dependency is not a completed review or accuracy pass.
The worker exits nonzero for a retained native failure. A successfully scored
cell still needs result readback/reporting and the next request's full
reviewed-prefix/fresh-capacity/held-release checks. QfO is explicitly refused:
it retains its separate native conversion, assessment and admission route.

All94 joined tests pass, zero failures/errors/skips,1.73s, including14 new
orchestration/CLI/envelope tests. Bash syntax and matching-runtime CLI import
pass; all920 frozen helpers remain unchanged. These are unit/preparation
checks, not successful real terminal execution. No new native score is claimed.
The original91-case receipt precedes the added source/CLI pin controls and is
preserved; the latest [JUnit](results/native_orthobench_followup_tests_20261005_v2.xml)
is the final tested worker version.

The dependent job's duration and resource use are postprocessing, not native
inference costs. Existing shared-host timing limitations stay intact: CPU,
memory-bandwidth and I/O competition has unknown potentially tool-dependent
effects. No quiet host, DGX, background subtraction or isolated ranking.

Keep the full publication goal active. Remaining native ablations,
uncertainty, biological strata, provenance and distribution requirements
are not completed by preparing this handoff.
