# Owned Preflight Worker Lifecycle

The local executor now starts and manages one environmental responder for its
single authorized run identity. The preceding
[responder implementation](THREADRIPPER_ENVIRONMENT_WORKER_20260928.md) supplied
the decision logic but required separate orchestration. This change connects
it to `run_threadripper_scaling.execute`; it does not submit a Slurm job or
create an approved policy/readiness record.

## Execution Order

Before creating an attempt session, the executor requires the readiness record
to contain an environmental-policy reference and checks that policy's schema,
review decision, host, frozen-plan identity and evidence. The responder still
checks the complete policy and its direct readiness binding independently.

After the existing runtime setup, the executor launches exactly one responder
as an ordinary child within its allocation. Its working directory is the
repository root; helper imports use that root rather than the native runtime's
PYTHONPATH. It disables user-site and bytecode writes, has a distinct unused
cache path, fixes PYTHONHASHSEED, and removes inherited Python-home and dynamic
loader overrides. These settings apply only to the responder, not to native
scientific commands. It inherits the job cgroup; no separate allocation is
created. The responder itself checks actual job membership.

The child starts before input preparation so it can wait for the guard's
parked-worker request. The guard uses the owner's file waiter: a child that
exits without a response fails promptly, and a response arriving near child
exit is rechecked before declaring it missing. Atomic response publication
and scientific/environmental validation remain the responder and guard's
responsibility; the owner does not treat any existing JSON as a pass.

After the guard binds the review, its budget callback first waits for a
successful responder exit. A nonzero exit or five-second join timeout prevents
the Slurm budget check. Only then does the existing fresh budget guard run,
followed by the collector's unchanged final observation/freshness check and
native release. Thus a still-running preflight responder is not deliberately
overlapped with inference. No collector or scientific-command settings changed.

## Failure And Evidence

The session retains `environment_worker.log`, an initial PID/command/environment
record, and `environment_worker_lifecycle.json` with terminal status, exit code,
parent exception type and any cleanup error. The executor's final result binds
the lifecycle file by hash. The policy reference is rechecked at release and
after measurement alongside the other frozen execution evidence.

On preparation failure, parent interruption or response timeout, cleanup signals
only this owned Popen child. It tries termination, waits five seconds, then uses
kill and a ten-second wait if necessary. It does not signal process groups,
enumerate targets to kill, cancel unrelated jobs or restart the responder.
If cleanup cannot establish termination, the record says so and the attempt
is not complete; an observation timeout does not prove the child stopped.
An existing parent exception is preserved while cleanup errors are recorded.
Normal return without a successful join is also an error, not a completed
preflight. An owner cannot be entered twice.

The child's preparation CPU/memory belong to the same job. Do not subtract
these costs or reinterpret the native-step memory peak as the whole-job peak.
Complete job accounting through final validation and teardown remains a
separate requirement.

## Validation And Limits

173 focused lifecycle, driver, responder, typed-process, policy, observer and
PSI tests pass. Lifecycle tests launch real, explicitly synthetic subprocesses
and verify response handling, successful join, nonzero exit, missing response,
response timeout, simulated join timeout, spawn failure, parent interruption,
termination and termination-resistant child cleanup. All test children are
terminal and reaped. No unrelated host process is signalled by these tests.

Driver composition tests verify start/measurement/join/budget/cleanup order,
that a failed join prevents the budget callback, retained failure evidence,
and refusal to reuse the attempt. These tests mock native measurement and the
worker in the driver; the real file-based guard/responder test and real child
lifecycle tests are complementary component checks, not a full native run.

No production identity, reviewed ordinary-service policy, numerical background
limit, readiness approval or execution permit was created. A complete native
handoff in the intended Slurm allocation, actual deadline behavior, whole-run
policy application, full-scale resource validation, final source recipe and
a verified quiet window still remain. The latest host-load evidence remains
the earlier busy-host capture; this turn did not perform another host poll.
The publication goal and all 27 planned timing identities are unchanged.
