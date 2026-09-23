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

This module has not yet been deployed or exercised against the DGX. The
previous short real-host shell probe is separate evidence, not validation
of this implementation. Integrate and test the held-session submission and
terminal-capture lifecycle, interruption recovery, periodic monitoring,
whole-host workload evidence, and frozen policy/recipe before scaling runs.
Abrupt session loss or SIGKILL cannot be repaired by a Python cleanup path;
the caller must recover the existing job identity and inspect actual service
and scheduler state without resubmission. Do not claim gap-free isolation.
