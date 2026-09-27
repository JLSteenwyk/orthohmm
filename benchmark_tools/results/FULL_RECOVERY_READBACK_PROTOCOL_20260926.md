# Full Recovery Readback Gate

The new `readback_full_recovery_orthobench.py` is bound to native job 22326,
the previously committed full-run plan and the prespecified protocol. It is
a new file and does not modify any source pinned by the active native run.
Thirty-eight focused readback, launcher and score-comparison tests pass.
The [live negative gate](full_recovery_readback_live_gate_20260926.json)
records an actual rejection of job 22326 while RUNNING, before scientific
readback or score generation.

## Admission Requirements

- Unique scheduler COMPLETED state with exit 0:0 and 32 allocated CPUs.
- Exact plan/protocol/submission/command and execution-start/completion chain.
- Exact isolated native interpreter, frozen scientific arguments and recorded
  environment, one attempt, no checkpoint reuse or production-default change.
- All native source, input and tool records covered by the original plan.
- Unchanged complete installation audits before and after native execution.
- Exact native and GNU time logs; no native failure receipt.
- Unchanged full input universe, all frozen reference files and scorer.

After admission, the existing four independent scientific readers check
structure, sequence content, reconciliation events/pairs and hierarchy. The
strict partition reader requires every input gene exactly once. The historical
score is recomputed and must match its retained complete score object before
comparison with the current result. Partition and all 70 family-score changes
are reported, not hidden or treated as a reason to tune parameters. Admission
is repeated after readback to detect changed artifacts.

```bash
python -m benchmark_tools.readback_full_recovery_orthobench \
  --repo /absolute/repo \
  --directory /absolute/repo/benchmarks/work/publication_full_recovery_orthobench_20260926 \
  --job 22326 --output /absolute/new/readback
```

The planned scheduler readback uses two CPUs, 64 GiB and a four-hour limit,
with an `afterok:22326` dependency and no requeue. It never launches or retries
native inference. A readback failure is retained for diagnosis; it does not
authorize a native rerun. Tests and queued readback are not full-data results.

## Queued Execution

After readback commit `548576c` was pushed, submitted
[job 22327](full_recovery_readback_submission_22327.json). Its readback plan
SHA256 is `cb90ce1254056bc60153b31dd1bdf506e4d695c2fd761175272f716f83cff78d`.
The bootstrap verifies this plan, the new readback source, native plan and
interpreter version before invocation, and checks the new source again after
completion. The original native plan independently pins existing harness
sources. Initial scheduler inspection confirms the audit PENDING on its native
dependency while native job 22326 remains RUNNING. This is not score admission.
