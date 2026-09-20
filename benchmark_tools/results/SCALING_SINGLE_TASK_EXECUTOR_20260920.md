# Replacement Scaling Single-Task Executor

Implemented `run_root_context_scaling.py` and its matching
`run_dgx_root_context_scaling.sh` entry point. The batch resource directives
match the allocation checker: exclusive20CPU/96GiB on spark-7ff0, one task,
24-hour limit, no requeue. The existing measurement adapter retains the
85800-second native limit and the original27-task plan unchanged.

## Launch Gates

The Python entry point rejects execution outside the intended deployed
recipe root or with a different recipe manifest path. It loads the pinned
v2 plan, checks all local Python helpers, the batch script and plan against
the recipe inventory, and reuses the native interpreter/host/affinity/cache/
loader preflight. The batch script additionally validates its four arguments
and rejects existing or dangling observer-cache paths before Python starts.

Execution requires a separately hash-pinned JSON authorization with:

- `schema`: `root_context_scaling_authorization_v1`;
- `execution_authorized`: boolean true;
- exact `plan_sha256`, `recipe_sha256`, integer `index` and integer `job_id`;
- `environment_policy` and `environment_preflight`: path/SHA-256 references.

The policy record must have status `approved_frozen_environment_policy` and
a nonempty `approval_reference`. The preflight record must have status
`environment_preflight_passed`, the same integer job identity, the policy
digest and boolean `whole_run_observer_ready=true`. All three documents are
checked for byte stability before and after measurement. These are external
gate receipts, **not a mechanism for granting approval or establishing the
truth of environmental claims**. Merely constructing matching JSON does
not approve a policy or make an observer operational.

No real authorization, approved policy, environmental-preflight producer or
whole-run observer deployment is provided by this change. The current
plan's `execution_authorized=false` is preserved; only a separate explicitly
approved per-job authorization can permit execution. The environmental
decision and user-service approval remain unresolved. Do not manufacture
passing receipts to launch. A future orchestrator must arrange the real
per-job preflight/authorization before invoking this entry point; that
orchestration has not been implemented or tested here.

## Execution And Evidence

After authorization, the executor reserves a fresh
`scaling_root_context_v1/sessions/task_XX` directory. A second invocation
for the same slot refuses to overwrite it. It queries `scontrol show job
JOB_ID --oneliner`, retains command/status/stdout/stderr and observation
timestamps, and requires the matching RUNNING allocation before calling
`measure_task` once. It does not submit, cancel, release, requeue or retry jobs.

The native timer remains inside the existing adapter, excluding controller
preflight and receipt writing. Preparation, runtime checking, copying,
native measurement and verification retain their existing stage boundaries.
The launch receipt pins gate documents and executor source; the result
records the wrapper status and verification-file identity. Exceptions and
interrupts save a failed receipt and propagate without retry. An uncatchable
termination can leave a launch without a result; that is incomplete evidence,
not authorization to restart.

The CLI returns zero only for a `command_exited_zero` wrapper status.
Its result still leaves scheduler-terminal verification, native-output
validation, environmental validity, timing admission and permission for
the next submission false. A failed native command is retained through the
wrapper status; raw replay and failure classification remain separate.
No source, service or output from an existing scientific run is changed.

## Verification And Remaining Integration

304 focused tests pass across the executor, allocation checker, measurement
adapter, task binding, raw measurement audit and controller recorder. Bash
syntax validation passes. New cases cover constructed authorization and
environmental receipts, denied/type-mismatched gates, source-byte drift,
fresh-directory behavior, controller mismatch, native exception/interrupt,
post-native gate mutation and CLI exit handling. These tests mock native
execution and controller responses; they are not a DGX end-to-end run.

The current development checkout is explicitly tested to fail the deployed
recipe gate. No recipe has been deployed and no replacement task submitted.
Remaining work includes approved environmental policy and genuine observation
collection, the bounded submitting session and terminal-record integration,
recipe freeze/deployment, and actual composed validation. Independent raw
replay, native outputs/failures and environmental evidence must be audited
before any controlled comparative resource claim. This executor alone does
not complete the publication objective.
