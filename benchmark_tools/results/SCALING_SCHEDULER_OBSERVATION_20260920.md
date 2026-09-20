# Long-Run Scheduler Observation

`capture_job_scheduler.py` previously exposed only its four-hour default
through the command line, although Python callers could set a longer bound.
The replacement scaling allocation can last 24 hours. The CLI now accepts
an explicit positive finite `--max-seconds`, while retaining the four-hour
default for existing callers. Example for a separately authorized run:

```sh
python -m benchmark_tools.capture_job_scheduler \
  --jobs JOB_ID --max-seconds 86580 --output NEW_CAPTURE_DIRECTORY
```

This 86,580-second observation example allows 30 seconds beyond the frozen
86,550-second local waiting-session limit; it does not change the native
timeout or scheduler allocation limit and does not authorize submission.
Queue delays consume the observation window. An absent terminal record is
retained as missing and must be investigated for the same job, not retried
as new native work. A live recorder must be confirmed before any future run
is released, and actual terminal/session/environment audits remain necessary.

Each controller query and inter-poll sleep uses at most the remaining
observation budget. The deadline is checked between jobs as well as passes.
The report records the requested duration, actual elapsed time and whether
observation stopped with jobs missing at its limit. OS scheduling, process
cleanup and final evidence writing can add latency beyond the requested
bound; this is not hard real-time execution. No cancellation or submission
command is issued by this recorder.

## Verification

The [live two-second smoke](scaling_scheduler_bound_smoke_20260920.json)
queried running local BLAST job 21713 once. It recorded zero observation
errors, `incomplete`, the same job as missing terminal evidence, and elapsed
2.00232 seconds. Process exit 1 is expected for incomplete capture. A separate
post-check confirmed the job still running. No scientific job was restarted,
cancelled or reconfigured. This is a local-controller smoke, not 24-hour DGX
validation or proof of future job lifecycle behavior.

Seventeen single-job recorder tests and twenty array-recorder tests pass.
The broader workflow suite passes all 1,139 tests in 46.40 seconds. New tests
cover invalid limits before filesystem changes, CLI forwarding, simulated
observation beyond four hours, remaining-budget query/sleep bounds, and
retention of nonterminal jobs without promotion or scheduler mutations.
No historical source recipe or archived measurement was changed. Frozen
environment policy, launch/session integration and scientific timing admission
remain incomplete.
