# Prospective Near-Full-Node CPU Controls

## Withdrawn Before Submission

This proposed experiment was redundant with the completed, independently
audited full-node panel 21918. That earlier panel and its descriptive analysis
were found during the implementation review; see
`FULL_NODE_CONTROL_RESULT_21918_20260919.md` and
`FULL_NODE_CONTROL_DESCRIPTION_20260919.md`. It already tested native steady
work, high process creation and known competing work in three balanced blocks.

No job was submitted under this proposal and no outcomes were collected.
The duplicate workload added at a017488 and its tests were removed; the
uncommitted duplicate collector was also discarded. The proposal below is
retained solely as a historical record, not an active prospective protocol
or justification for selective repeat experiments. Continue from the existing
audited results, retaining all native-tool flags and failed overhead tasks.

## Status And Question

Protocol prepared before collection. The bounded workload generator is
implemented and locally unit tested; the integrated DGX collector/launcher
has not yet been implemented or submitted. No DGX outcome is available.

The earlier diagnostic found 21 satellite_v2 narrow CPU flags inside its
phylogeny stage, associated with higher host process creation but not higher
named outside-job CPU. This association does not identify a cause. Test
whether heavy known native process creation can reproduce flags, while
known outside-native CPU remains detectable. These controls cannot recover
missing observations or admit the earlier timing/scaling panels.

## Fixed Design

Use a fresh exclusive DGX spark-7ff0 allocation with 20 CPUs and 96 GiB,
30-minute maximum, no requeue. No unrelated job/service may be stopped.
Check the scheduler before launch, then avoid interactive SSH while the
panel is running. Record controller observations separately. Preserve every
trial and failure; no selective repetition.

Use three complete blocks in this fixed balanced order:

| Block | First | Second | Third |
| --- | --- | --- | --- |
| 0 | steady | creation | positive |
| 1 | creation | positive | steady |
| 2 | positive | steady | creation |

For every trial, use a distinct native step and retain ready/go/done/release
handshakes. Require at least 20 allowed CPUs in the allocation and record
their identities/topology. Pin 18 native worker processes to distinct CPUs;
reserve a nineteenth for the observer and a twentieth for the positive
control. Require recorded memberships and affinities to match these roles.
Do not interpret 18 workers as all 20 cores or as equal compute capacity
across the heterogeneous DGX CPU types. Keep the same assignments in all
trials. Native workers inherit the native step's cgroup; positive control
inherits the observer's batch cgroup, not the native step.

- **steady:** 18 workers perform bounded arithmetic for 20 seconds each.
- **creation:** 18 workers repeatedly fork one child at a time; each child
  targets 0.005 process CPU seconds of the same arithmetic before exit.
  Each worker runs for at most 20 seconds, with at most 10,000 children.
  Each child has a one-second monotonic work bound and cannot extend past
  its worker's work deadline. Wait/reap each child before creating another.
- **positive:** same 18 steady native workers, plus one 20-second steady
  worker on the reserved batch CPU, started after the initial observation.
  This is known outside-native demand, not an unrelated user's job and
  not necessarily contention on a native CPU.

`cpu_creation_workload.py` records per-worker monotonic start/end, affinity,
cgroup membership, self CPU, reaped-child CPU, child count, wall-limited
children and cap attainment. Its launcher bounds completion at requested
duration plus 15 seconds and terminates only its own worker process groups
on cleanup, including live fork children. Launch startup is not included
in the recorded worker CPU sum; report this rather than subtracting it.

## Collection And Interpretation

Use the existing dual-bracket counter reader and unchanged screens at
one-second cadence, with observations before work and after all workers
finish while the native wrapper remains alive. Retain both outer and narrow
screens, native PSI, frontier snapshots and all invalid-point evidence.
Any transient-cgroup invalidation stays a failed measurement, not a skipped
point. Freeze launcher, helper identities and this protocol at a source
commit before submitting; do not modify previously frozen collectors.

Require worker exit success, no child-count cap, matching role memberships,
and at least 15 seconds of simultaneous overlap among native workers.
For a near-full-node interpretation require at least 240 aggregate worker
CPU seconds per trial. For a creation witness additionally require at least
1,800 reaped children. For the positive witness require at least 12 CPU
seconds and 15 seconds overlap of the batch worker with the native common
window. If any witness fails, report the achieved dose and call the
corresponding diagnostic inconclusive, not a negative finding.

Report all intervals. Separately summarize intervals whose full observation
windows lie inside the common worker window; require at least ten such
intervals for a within-work interpretation. For the positive condition use
the common window that also includes its batch worker. This prespecified
interior summary must not replace or erase startup/teardown observations.

Primary diagnostic summaries per trial are the number/fraction of narrow
`excess_unassigned_cpu` flags and the median/maximum signed residual, with
the analogous outer-screen values. Report actual CPU, creation counts/rates,
PSI and read-window durations alongside them. Do not treat adjacent intervals
as independent replicates. Present three blockwise contrasts descriptively;
do not generate a p-value from interval counts or tune a threshold.

Native-only flags with validated creation activity would show that known
native work can coexist with flags, not that every historical flag is a
false positive. No creation flags would not rule out exec-heavy tool-specific
behavior. Positive detection requires at least 80% of eligible interior
intervals to flag excess unassigned CPU in each valid positive trial.
Failure of that check prevents a detector-sensitivity claim. Negative
residual discrepancies or steal increments remain separate reported flags.

These controls do not validate observer overhead, handle transient cgroup
coverage, exclude memory/I/O interference, prove host quietness, or admit
scientific timing. No wall-time correction, historical reclassification,
threshold relaxation or direct transition to the 27-run scaling panel is
authorized by this experiment.
