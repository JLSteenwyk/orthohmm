# Threadripper Collector Diagnostics

These are engineering observations, not scientific timing results. The
[receipt](threadripper_collector_diagnostics_20260928.json) pins both attempts'
raw files and the new collector/replay sources. Historical DGX collectors
and their replay rules were not edited.

Both jobs requested bizon, partition gpu, one exclusive node, one task,
64 scheduler CPU slots, 128 GiB RAM, three minutes, and no requeue. The
batch observer invoked `measure(["/bin/sleep", "35"], directory, job_id,
32, 128*1024**3, 85800, 1.)` using the base Python. The collector launched
a separate srun step with 64 slots and `--cpu-bind=mask_cpu:0xffffffff`.
Output directories are `benchmarks/work/threadripper_collector_probe_JOBID`.
No GPU was requested, no DGX was accessed, and no unrelated job was stopped.

## Outcomes

- **22341: FAILED 1:0, 39 seconds.** The native sleep completed, but final
  evaluation rejected the initial interval. The first two point starts were
  2.142260239 seconds apart because whole-host inventory ran between the
  first sample and the cadence start. All 37 raw points and the traceback
  remain retained. No complete collector report or timing admission exists.
- **22342: COMPLETED 0:0, 40 seconds.** After explicitly moving the initial
  host inventory before the first point, the unchanged interval evaluator
  accepted all 37 points. Replay reproduced the native duration
  35.031039147 seconds, root/lineage reports, memory evidence and host
  process evidence. All affinity observations reported CPUs within 0-31;
  these observations included the wrapper and its sleep child, not a
  multithreaded scientific engine.

The second attempt is a documented engineering correction, not an automatic
rerun or selection among scientific timing repeats. The 27 planned inference
identities remain unrun. No sampling thresholds were relaxed.

Three successful host-process samples bracketed the command with no read
errors. Maximum observed outside-job usage was **70.7703828113 CPU-core
equivalents**. This establishes observed competition, not a quiet window.
No background-work subtraction was applied.

## Verification and Remaining Work

126 focused tests passed in 37.02 seconds, covering the new collector settings,
scope selection, affinity races/violations and replay tampering, plus the
historical scaling collector/replay and root-context tests.

The new replay preserves observed affinity violations and gaps instead of
converting them into compliance. It does not authenticate acquisition by
itself or authorize scientific admission. Placement is inherited affinity,
not a hard 32-CPU quota; swap remains unlimited under the retained protocol.

Still required: live multithreaded and deliberate-widening controls, native
x86 runtime/input-enumeration validation, observer-overhead assessment,
executor/provenance integration and a verified quiet window. The collector
still retains all points in memory, so long-run observer resource use must
be assessed before production. This diagnostic does not resolve broader
scientific uncertainty, source retrieval or publication-release requirements.
