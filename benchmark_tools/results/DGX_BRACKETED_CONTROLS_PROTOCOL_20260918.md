# Bracketed CPU Positive Controls

Freeze this protocol and code before execution. One exclusive spark-7ff0 job,
two CPUs/task,2GiB,five minutes,no requeue. Execute three trials in fixed order:
quiet, completed one-process-CPU-second sibling burst, and sustained
three-process-CPU-second sibling load. In every trial a separate one-CPU
native Slurm step performs two process CPU seconds of integer arithmetic.
Do not run scientific inference or stop unrelated processes.

The native worker announces readiness without taking its pre-work snapshot.
The parent completes its host pre-snapshot, then releases the worker. The
worker takes its native pre-snapshot, performs fixed-CPU-duration work, and
takes its post-snapshot. The parent takes its post-snapshot only after both
the native step and any sibling load have exited successfully. The sibling
must remain in the batch scope, separate from the native step. Its CPU duration
must be in[requested,requested+0.1)seconds; native duration in[2,2.1)seconds.

Apply the frozen screen_bracketed_cpu thresholds: residual at most0.25cores,
negative tolerance0.5CPU seconds, no host steal. Expect quiet to pass and
both load trials to fail the operational screen. Retain every result even
if those expectations fail; do not change thresholds or retry selectively.
Retain raw snapshots, commands, source hashes, measured CPU durations, control
expectations and terminal scheduler evidence. A successful process exit means
the controls were evaluated, not that all expectations necessarily held.

The sustained sibling may outlast native work; the residual intentionally
includes this outer-window excess. It is not a correction to native runtime
or an estimate of CPU interference during native work alone. This is a
coarse positive-control experiment, not general calibration or statistical
uncertainty. It does not grant scientific timing admission, address non-CPU
interference, validate interval-level monitoring, or authorize replacement
scaling runs. Keep all publication/controlled-timing flags false.
