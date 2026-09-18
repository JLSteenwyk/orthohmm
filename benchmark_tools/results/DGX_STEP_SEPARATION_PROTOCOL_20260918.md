# Prospective Native-Step Separation Probe

Engineering only. Run one exclusive Slurm job on spark-7ff0,2CPUs/task,2GiB,
5-minute limit,no-requeue. No native OrthoHMM/OrthoFinder inference is run.
The batch observer starts a one-CPU srun step containing a sleeping worker.
File handshakes ensure the native step's first snapshot precedes each observer
window and its final snapshot follows it. Job IDs must match, while batch and
native step scopes must differ. Do not move or signal unrelated processes.

Execute exactly two trials in order: quiet (one-second sleep), then burst
(a child in the batch cgroup consuming at least0.75 process CPU seconds).
The burst must exit successfully before the observer's final snapshot and
before release of the native worker. Its reported CPU time must be below5s.
Require host accounted busy CPU increment of at least0.5s during the burst
window. This is a coarse positive-control check, not a calibrated discrepancy
bound; the quiet result is descriptive and has no admission threshold.

Retain both trials, raw native/batch snapshots, burst CPU/time/scope evidence,
commands, source hashes and terminal scheduler status. Preserve failures;
do not retry selectively or interpret an engineering failure as a method loss.
Write no complete success report unless both steps exit successfully and all
checks pass. The probe must leave publication_ready and
controlled_workload_verified false even on success.

This tests whether host cumulative accounting retains a known exited child
outside a distinct native step. The child remains within our allocation; it
does not test an unrelated job, general transient detection, observer overhead,
native compute/process-heavy behavior, memory peaks or scientific timing.
It cannot upgrade the original27 descriptive timing runs. A later timing
experiment still requires a prospectively frozen inclusion/execution plan.
