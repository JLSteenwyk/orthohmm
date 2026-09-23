# DGX Service Guard Integration Contract

`benchmark_tools.dgx_service_guard.ServiceGuard` implements the approved
Samwise stop, runtime mask, repeated state checks, and restoration policy.
It is a session component, not a submitter, watchdog, or timing authorization.
No real timing allocation was launched by this change.

The caller must keep one SSH session alive throughout the allocation. Call
`begin()` before submission, `before_submission()` immediately before the
single submission attempt, and `bind_job(job_id)` as soon as its identity is
known. Call `check()` periodically throughout the job. A lost mask, persistent
configuration drift, or failed state observation invalidates suppression
evidence and must prevent further benchmark submissions. Do not construct a
fresh guard to bypass an unresolved submission or an existing mask.

`restore()` queries the bound job with `sacct` and requires its exact raw ID,
a terminal state, and a parseable terminal timestamp. Failed and timed-out
jobs also require restoration; successful execution is not a prerequisite.
An empty result, running/completing state, query failure, timeout, or unknown
submission identity leaves suppression in place. A polling timeout is never
interpreted as job termination. Before any submission attempt, prelaunch
errors can instead be cleaned up without a scheduler query.

Only the exact approved service is controlled. Persistent unit bytes and its
enabled link are checked before changes and restoration. Existing runtime
masks are rejected, and only the owned `/dev/null` mask is removed. Originally
inactive services remain unstarted; otherwise prior start behavior is
requested. This is not an application-health claim: scientific_openclaw is
still missing. Command receipts and restoration evidence use fresh paths.

## Verification And Remaining Work

The combined guard, scaling-executor, and existing session-submitter unit
suites passed **81 tests**. Guard cases include each accepted prior state,
terminal failure states, live/missing/wrong/duplicate accounting, malformed
timestamps, query timeouts, submission uncertainty, lost/changed masks,
configuration drift, and prelaunch stop-observation failure. These use a
fake systemctl/sacct transport and real temporary filesystem paths.

The initial implementation had not been deployed; the real-host follow-up
below now validates pre-submission service control only. Integrate and test
the held-session submission and
terminal-capture lifecycle, interruption recovery, periodic monitoring,
whole-host workload evidence, and frozen policy/recipe before scaling runs.
Abrupt session loss or SIGKILL cannot be repaired by a Python cleanup path;
the caller must recover the existing job identity and inspect actual service
and scheduler state without resubmission. Do not claim gap-free isolation.

## Real-Host Prelaunch Checks

Deployed the exact guard, `probe_dgx_service_guard.py`, and its two helper
modules under `/home/jlsteenwyk/projects/orthohmm-publication/service_guard_probe_20260923`.
The spark queue was empty before execution. Each invocation used one held
SSH session, explicit approved-service environment marker, system Python
3.12, and a fresh receipt directory. No Slurm job was submitted.

- Normal probe: three suppression samples at two-second intervals; 11
  successful systemctl receipts; restoration verified.
- Deliberate SIGTERM after the second sample: exception cleanup restored
  prior start behavior; 10 successful systemctl receipts.
- Independent local receipt inspection confirmed all command exits zero,
  every intermediate service state masked/inactive/PID0, and final loaded
  state. Both prior states were activating. Copied and current source hashes
  matched every recorded source identity.
- Subsequent SSH inspection confirmed loaded/activating/auto-restart with
  ExecMainStatus 1. This is restoration of the original broken application's
  behavior, not successful installation or health verification.

Preserved command, sample, prior-state, restoration and probe receipts in
`dgx_service_guard_receipts_20260923.tar.gz`, SHA-256
`dbd268367878a2b4ee29f55426d10e8e3a4023fc43eb9cd330292471fa487a1f`.
The guard source SHA-256 is
`a67a6e6f90327da6841d8c62b3248f5419447d7b3498956158c745c95f00a32c`;
the probe source SHA-256 is
`3370dc5c55ad42d893e8040c314052a435a85179f274b71cf04af25954fadf03`.
The guard's 34 unit cases passed again before deployment. These short
prelaunch tests do not validate bound-job terminal gating on the DGX,
arbitrary SSH loss, long-duration masking, or whole-host isolation.
