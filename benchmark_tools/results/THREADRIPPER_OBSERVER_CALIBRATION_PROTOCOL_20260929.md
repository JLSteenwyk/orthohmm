# Threadripper Descendant Accounting Calibration

## Prospective Scope

This is one bounded engineering diagnostic on `bizon`, not one of the 27
production scaling identities. Run after job 22377 becomes terminal, regardless
of its outcome (`afterany:22377`); the diagnostic does not consume its results.
Do not stop unrelated work. Shared-host contention is retained and prevents
interpreting this run as controlled overhead or comparative efficiency evidence.

The [machine-readable protocol](threadripper_observer_calibration_protocol_20260929.json)
is frozen before submission at SHA-256
`ebe0aabcaf32111a6035aacd7a15d03a3197f47fca9a9f722ee0ad5be5c8fb7e`.
It pins 822 Python harness sources, the private controller interpreter, the
submission script and the already-adopted resource endpoint protocol. No frozen
scientific inference source or endpoint helper is changed. Source changes before
launch fail validation; do not silently refresh pins or retry the same attempt.

## Workload And Checks

- 32 independently launched Python worker processes, four additional threads
  per worker, 64 MiB retained/page-touched allocation per worker, 30 seconds
  of busy work. All inherit CPU affinity 0-31. Python threads share a GIL:
  this is process diversity and thread enumeration, not 128-core scaling.
- Slurm: one exclusive `bizon` allocation, 64 CPUs per task, 128 GiB RAM,
  26-hour allocation ceiling compatible with the unchanged collector timeout;
  no requeue. The actual diagnostic workload is bounded to 30 seconds plus
  readiness, collection and cleanup. Fresh output only; one attempt.
- Require successful native exit and at least five complete interior samples
  containing all 160 witnessed worker/main thread identities, matching start
  ticks, cgroups and affinity. Every sample wholly inside the shared worker
  interval must be complete; races or missing threads fail this check.
- Require each raw point read to take at most one second, and interior sample
  starts to be at most 1.5 seconds apart. These are cadence feasibility checks,
  not a numerical bound on causal workload slowdown. File writing, host scans,
  reporting and scheduling costs are not all represented by point read time.
- Sum worker user/system `getrusage` interval deltas independently. Native CPU
  bracket accounting must contain that sum within 0.01 seconds of rounding,
  with excess at most the greater of 10 seconds or 5% of witnessed worker CPU.
  This permits uncorrected startup/wrapper work; it is not an overhead
  subtraction or a precision claim about pure algorithm CPU.
- Native-step lifetime kernel peak must contain the simultaneously retained
  2 GiB allocation. This is a lower-bound check, not exact RSS calibration;
  launcher, cache and other native-step memory remain included. Whole-job
  reporting peaks stay supplementary and overlapping peaks are not added.

The worker readiness gate has one global 35-second deadline. Failed startup
terminates and reaps only owned children; partial logs remain. Workers retain
PID/thread start ticks, memberships, affinities, allocations, monotonic times
and process CPU witnesses. No failed or incomplete attempt is automatically
replaced. A failure requires diagnosis and a separately justified new protocol.

## Independent Readback

After scheduler completion, inspect allocation/batch/native-step exit statuses
and retained logs. The independent audit replays raw collector evidence and
re-derives the adopted scoped resources rather than trusting the diagnostic's
precomputed summary:

```bash
python -B -m benchmark_tools.audit_threadripper_observer \
  --directory benchmarks/work/threadripper_observer_calibration_20260929 \
  --protocol benchmark_tools/results/threadripper_observer_calibration_protocol_20260929.json \
  --protocol-sha256 ebe0aabcaf32111a6035aacd7a15d03a3197f47fca9a9f722ee0ad5be5c8fb7e \
  --output benchmarks/work/threadripper_observer_calibration_20260929/independent_audit.json
```

Use a fresh output. Failed numerical checks are saved with exit status 1;
invalid provenance/input structure fails rather than granting acceptance.
Retain failures and partial derived receipts. Passing checks do not establish
continuous containment, a quiet host, native full-workflow handoff, causal
observer slowdown, scientific timing admission or publication readiness.
Those remain separate requirements before production timing. No DGX operation
or scheduler/service configuration change is authorized by this protocol.

## Validation Before Submission

60 focused tests passed across the new diagnostic/auditor, scoped endpoint
derivation and raw replay. Coverage includes bounded dimensions, malformed
worker indices, a real two-process/two-thread/1-MiB/one-second smoke run,
refusal of repeat output, owned-child cleanup on startup failure, missing or
reused thread identities, scope/affinity races, cadence, CPU and memory checks,
and protocol source/envelope validation. The smoke run is not the scheduled
full diagnostic and supplies no comparative timing evidence.
