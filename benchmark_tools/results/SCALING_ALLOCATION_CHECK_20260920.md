# Replacement Scaling Allocation Check

Added `benchmark_tools.verify_scaling_allocation` for the per-job executor
and subsequent audit. Environment variables alone do not prove controller
allocation settings. This read-only checker accepts exactly one oneline
controller record for the specified non-array, non-heterogeneous job and
checks the replacement panel's 20 CPU, 96 GiB, single-node/task, exclusive,
non-requeued, zero-restart, 24-hour allocation on spark-7ff0/partition spark.

The expected entry point is
`root_context_scaling_recipe_v1/benchmark_tools/run_dgx_root_context_scaling.sh`
under the frozen DGX project root, with that recipe directory as WorkDir.
This names the prospective executor contract; the script and deployment are
not supplied or claimed complete by this change. No recipe or frozen plan
was modified. A future executor must use this path or explicitly revise the
prospective contract before execution.

`running` accepts only RUNNING, not PENDING or COMPLETING. `terminal` requires
a terminal controller state and retains its exact exit/signal fields,
including failures, rather than equating termination with native success.
An expired observation cannot become terminal evidence. The two equivalent
24-hour printed limits are accepted; shorter/longer/unlimited allocations
are rejected. Duplicate fields, multiple records and missing required fields
are rejected rather than silently selected.

The file audit binds the pinned v2 plan and scheduler bytes and rechecks
both after validation. A requested index selects a valid frozen slot but
**does not bind that slot to the job**: that still requires the submitted
arguments, recipe, task records and session receipts. The result explicitly
retains `task_identity_bound=false` and all execution/continuation/timing
authorization flags as false.

Read-only invocation from the repository root, using an already retained
controller record and a fresh report destination:

```bash
python -m benchmark_tools.verify_scaling_allocation \
  --plan benchmark_tools/results/dgx_root_context_scaling_plan_v2_20260920.json \
  --index 0 --job JOB_ID --phase terminal \
  --scheduler /path/to/retained/scheduler_JOB_ID.txt \
  --output /path/to/new/allocation_audit.json
```

`JOB_ID` and paths are placeholders, not a submitted replacement run. This
command does not query, submit, cancel or retry a job, alter a service, or
run inference. The pure validation API can also check a freshly captured
RUNNING record before inference. Freshness/authenticity of that capture is
the caller's responsibility, not inferred from the record text.

Validation: 199 tests passed across this checker, frozen task binding,
single-job controller capture and composed measurement replay. Constructed
controller fixtures cover every retained terminal state, all 27 frozen slot
indices, policy mismatches, missing fields, type/phase errors, mutated input
evidence and CLI no-overwrite behavior. These are not 27 executed native
runs or a real accepted DGX24h allocation. No replacement job was launched.

The environmental-policy freeze, authorized executor, bounded session and
terminal capture integration, recipe deployment and actual end-to-end
validation remain required. Exclusive Slurm allocation does not exclude
user services or other work outside Slurm and does not admit comparative
resource results. The publication objective remains incomplete.
