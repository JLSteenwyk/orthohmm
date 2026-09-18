# Complete-Command Interval Native Smokes

## Frozen Scope

Run the same three native commands and645-protein/eight-species fixture as
the retained counter-native smoke, with the original smoke specification
SHA-256`8af8a049c69b7e479e8443e1e0e14b74c78d67b345c03d7f557768488d057fa2`.
Relocate only the output prefix to`interval_native_smoke_v1`. Retain scientific
flags, seeds, environments, original inputs and native enumeration. No
scientific parameter changes or timing comparison with old smokes are planned.

Use a sequential%1 array of three exclusive spark-7ff0 tasks,20CPUs/96GiB,
one-hour allocation limit and900-second native timeout, no requeue. Freeze
and push implementation/protocol and verify a fresh remote recipe manifest
before execution. Recheck runtime trees, recipe and inputs before/after each
native measurement through the existing verified launcher. These checks and
input preparation remain outside the command timer.

## Measurement

Worker publishes PID/cgroup before its initial counter snapshot and waits.
Observer records a complete host/native-step/host observation point, then
releases the worker. Worker records native-before counters, executes the
unchanged GNU-time-wrapped command, records native-after counters and publishes
completion. Worker then waits for release, preserving its step counters.

Observer takes points on a one-second schedule. It recognizes completion
before taking the final point, guaranteeing that the point encloses the
complete command. A completion arriving during a point is recognized on the
next scheduled point; never invent a shorter final interval. Reject missing,
overlapping or irregular intervals with the frozen interval validator.
Before releasing the worker, read step-aggregate memory.current, memory.peak
and memory.events with scope/read timestamps and retained read errors.

Retain native command exit and timeout status. On timeout kill only the
command's own process group; Slurm also bounds the allocation. Failed commands
and observations remain in their output directories and are not selectively
rerun. Whole-command screening uses native-before/after snapshots and enclosing
outer host reads. Interval screening uses the native step aggregate. Preserve
both scopes explicitly. Neither interval residuals nor native wall time are
corrected or summed. Diagnostic flags do not suppress valid native outputs.

## Validation And Limits

Require full boundary/interval replay, before/after runtime and input checks,
completed scheduler status, unchanged relocated commands, native partition,
pair and graph validation, and finite GNU-time output. Retain screen flags,
memory-read failures and other adverse observations rather than presenting
only quiet runs. The step memory peak includes wrapper/startup and cache and
is not process RSS.

One fixture per method checks integration, not accuracy, scaling or general
overhead. Observer and native scopes share allocation CPUs. Accounting delay,
non-CPU interference, CPU frequency/thermal state and a prospective scientific
inclusion/repeat policy remain unresolved. No original timing is upgraded and
no27-run replacement panel is authorized by this experiment. All scientific
timing and publication-ready flags remain false even if all checks pass.
