# Single-Identity Threadripper Executor

`benchmark_tools/run_threadripper_scaling.py` joins the pinned local command
plan, lookup/runtime checks, reviewed panel history, input preparation,
collector and scheduler release-budget guard. The companion submission script
requests the unchanged exclusive 64-slot task, 128 GiB and 26-hour allocation.
Native commands retain 32-CPU affinity and the 85,800-second timeout.

This is an execution component, **not production readiness or environmental
certification**. No real readiness review, execution permit or production
submission was created. Diagnostic fixture success and an empty queue cannot
satisfy the missing prerequisites.

## Required Inputs

The CLI takes `--request ABSOLUTE_JSON --request-sha256 SHA256`. The immutable
request must bind the actual scheduler job ID, next index, local allocation
cwd and submission script, exact plan/lookup digests, executor recipe,
ordered prior session references, readiness review and run-specific future
environmental-review path. Requests for skipped, repeated, live, unreviewed
or otherwise ineligible identities are rejected by the existing history
binder and panel state machine. It does not independently authenticate the
historical review conclusions; that remains the external audit's task.

The executor recipe inventories every top-level `benchmark_tools/*.py` file
and the submission script, with direct source records. It must be frozen
after source stabilization. No changed helper can be silently omitted.
The readiness review binds that recipe and the current plan/lookup, cites
supporting evidence and an explicit review reference, and requires actual
full-scale observer validation and a frozen environment policy. Boolean
fields are assertions to verify, not substitutes for these analyses.

Job-specific authorization must be constructed after obtaining a real job ID
and before allowing the submission to execute, for example using a held
allocation. The script does not submit, release, cancel or retry any job.
It isolates interpreter bytecode lookup at startup; no existing cache is
deleted. The executor creates one fresh `sessions/run_XX` directory and the
existing wrapper creates the frozen native run directory. Existing attempts
are never overwritten. Pre-entry failures remain scheduler/log failures and
must also be retained by the external submission ledger.

## Environmental Handoff

After preparation and pre-native runtime checks, the collector parks its
native worker. The release guard writes `environment_review_requested.json`
in the measurement directory. A separate review worker must inspect fresh
host/scheduler evidence under the already frozen policy, then atomically
publish the requested session's `environment_preflight.json` within 20 seconds.
This worker is not implemented or approved by this change. A subsequent
[process-policy comparison component](THREADRIPPER_PROCESS_POLICY_20260928.md)
has unit and retained-fixture negative checks; it does not establish a reviewed
service inventory or complete the live worker and handoff.

The preflight must identify this job, index, recipe and readiness review,
include supporting file records and a review reference, and report the
whole-run observer ready and no unrelated scientific workload. Its observation
must start after this release request, belong to the same boot and be no more
than 120 seconds old. A preexisting review, changed evidence, missing response,
future timestamp or negative decision prevents native release. The accepted
review's full file identity is recorded at handoff; it is not retrospectively
inserted into the earlier request. The source/history records are rechecked
before the request. The live scheduler budget check follows this review.

The 20-second response bound leaves room within the native worker's unchanged
45-second gate timeout; it does not guarantee success on an overloaded host.
Timeouts fail closed and remain failed attempts. No retry or automatic next
submission follows. Binding an external review does not establish that its
interpretation of the host is correct. The existing monitor's 0.25-core
diagnostic threshold has not been converted into an eligibility threshold.

## Validation And Remaining Work

141 focused tests pass across the executor, history binder, panel state
machine, controller guard and local measurement wrapper. Tests include source
drift, wrong jobs/indices, incomplete readiness, stale/wrong-boot reviews,
missing/preexisting handoffs, timeout, failure retention, refusal to repeat
an attempt and a real subprocess/file handshake with explicitly synthetic
review data. Shell syntax validation passes. The composition tests mock native
measurement; they are not integrated native or production timing validation.

The [fresh host observation](threadripper_executor_host_observation_20260928.json)
finds about 104 observed competing CPU-core equivalents despite an empty
Slurm queue, including IQ-TREE, Python, BAli-Phy and gainLoss work. No unrelated
process was signalled. A quiet window is not established.

Before production: finish full-scale observer/resource validation, freeze and
review the ordinary-service policy, implement/test its live review worker,
freeze this executor's recipe, and validate the complete handoff with the
native collector. Do not manufacture passing review files to bypass these
requirements. After every actual run, independently audit terminal scheduler
state, runtime, environment, resources, and native outputs or retained failure
before creating its reviewed history session and considering the next index.
The original 27 identities/order and scientific settings remain unchanged.
