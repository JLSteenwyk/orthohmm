# Prospective Incremental Root-Context Collector Comparison

## Question And Frozen Work

Estimate the observed incremental elapsed-time change from adding the
separately timestamped root/system/user/init context reads to the existing
periodic lineage collector. This is a paired engineering comparison, not
an accuracy benchmark or a method-speed ranking. Native integration 22021
passed, but it did not measure this incremental overhead.

Derive every native command from the pinned three-task root-context native
plan (SHA-256 `9b99cd810aaf7e040cda0242dd9b6d1da82bcb231c22657b8dc0108c4c4178ad`).
Preserve all 73,266 proteins, four proteomes, method arguments, source/runtime
identity, enumeration, environment, worker allocation and GNU-time wrapper.
Only fresh output/cache paths, collector selection and diagnostic metadata
change. Do not use historical runs as the paired baseline.

## Fixed Design

Use one exclusive DGX20CPU/96GiB allocation, serial execution, no requeue,
and a five-hour job limit. Each native step retains a 900-second timeout.
There are 18 tasks, nine adjacent pairs, with three pairs per method:

| Block | Method order | Order within every pair |
| --- | --- | --- |
| 0 | high sensitivity, satellite_v2, full OrthoFinder | lineage, root context |
| 1 | satellite_v2, full OrthoFinder, high sensitivity | root context, lineage |
| 2 | full OrthoFinder, high sensitivity, satellite_v2 | lineage, root context |

The lineage arm uses `measure_native_lineage_step.measure`; the root-context
arm uses `measure_native_root_context.measure_native_run`. Both are periodic
at one second, with the same native worker/timer and host-monitor settings.
This is not a boundary-only-versus-periodic comparison. Both arms retain
original/narrow lineage screens; the added root-context report is present
only in the root-context arm. Keep the method configuration fixed within
each pair and across blocks.

Freeze the plan and complete source recipe before submission. Require an
empty DGX queue immediately before the single submission. Hold the submitting
SSH session until termination with remote/local bounds of 18120/18150 seconds
and a 10-second remote termination grace. Record bounds before launch. Do not
open additional DGX SSH sessions or transfer files during the panel. Preserve
unrelated services and login settings; retain the bounded user-manager journal
afterward. Observation timeout is not terminal evidence and must not trigger
an automatic resubmission.

Stop subsequent tasks on preparation, native execution, runtime verification,
measurement or cleanup failure. Retain the failed task and all unrun entries.
Do not retry failed pairs, change method settings, suppress flags, select a
clean subset, or infer missing outcomes. A post-run output mismatch remains
invalid paired evidence even if the executable exited successfully.

## Endpoints And Decision Rules

Audit all native output semantics and exact input/runtime/source/command
bindings. Replay the appropriate raw collector in each arm. Require canonical
output equality within each pair and against the validated 22021 same-method
outputs, whose audit SHA-256 is
`9b2a243a13cd83ab431025a7f24d5b166672d50cbaa581697a2eba2557e54b23`.
Check within/across-task observation ordering and scope identity, and retain
GNU-time CPU, native-step cgroup memory/pressure and all flags.

For each valid pair, compute the signed ratio
`root_context_native_wall_s / lineage_native_wall_s - 1`. Do not clamp negative
values or subtract this estimate from any scientific runtime. Report all nine
pair differences and the median of exactly three valid pairs for each method.
If a method lacks any valid pair, its three-pair median and budget conclusion
are null, not recomputed from the remaining subset. The overall panel budget
conclusion is null unless all nine pairs and panel provenance are valid.

Retain the prior engineering budgets: every valid pair's signed ratio at most
0.10 and each complete method's median at most 0.05. These are engineering
thresholds, not confidence limits, a population overhead bound, or a claim that
negative observations imply the collector accelerates computation. Report
order/block effects descriptively; three pairs are not a precision estimate.

## Limits And Next Gate

The order is fixed, not randomized, and collector order is not perfectly
balanced within a method. Time drift, cache effects and background-service
changes may affect paired differences. The waiting session and unrelated
services are part of the environment. Passing the elapsed-time budgets does
not remove flags, establish causal attribution, exclude non-CPU interference,
or admit scientific scaling timings. Failure remains a result, not grounds
for endpoint changes. No method defaults or accuracy scores change here.

This experiment concerns incremental root-context reads on the existing
four-proteome workload, not total instrumentation overhead or larger-input
scaling. Before scientific timing admission, explicitly assess the remaining
environmental and measurement limitations against the publication protocol.
