# Typed Whole-Run Process Evidence

The local Threadripper collector now passes `enriched_snapshot` to the existing
streaming `HostMonitor`. Initial, periodic and final process records therefore
retain same-boot Linux Kthread evidence and identity errors, not just CPU
counters and names. The historical shared host monitor and DGX collectors are
unchanged. The native command, input bytes, CPU allocation, 30-second target
host cadence and 27-run scientific order are unchanged.

The additive snapshot fields retain the v4 resource report format. They do not
retroactively upgrade old observations. Helper-source recipe hashes change and
must be frozen again before native validation or production. Typed reads add
observer work and latency; previous native fixtures do not validate that cost.

`benchmark_tools.review_threadripper_process_stream.evaluate` replays the raw
JSONL stream with only the previous and current snapshot retained in memory.
It ignores stored interval verdicts and recomputes every adjacent comparison
against the explicitly supplied v2 process policy. It requires consecutive
indices, a consistent observer, native-interval bracketing, and caller-supplied
prospective limits for foreign CPU and start-to-start sampling periods.
Snapshot duration also must fit the period bound. A missing, malformed or failed
record breaks the chain; later good observations cannot erase the failure.
Every interval is checked, including intervals surrounding native execution.

No numerical policy defaults were selected. The new period bound is a future
policy requirement, not an inference from an observed runtime. The function
does not bind file hashes itself: callers must verify raw input, policy and
native-interval provenance before using its result. Its positive status is only
`sampled_process_policy_satisfied`; both controlled-workload verification and
scientific timing admission remain false.

## Validation And Limits

204 focused process-stream, collector, process-policy, typed-identity, host
monitor, executor, environmental-worker, replay, wrapper and child-lifecycle
tests pass. The tests include real HostMonitor serialization with synthetic
typed observations, a one-pass 1,000-snapshot stream, and rejection of
middle-only competitors, CPU bursts, changed identity/type/boot, missing
records, wrong observer, counter decreases, oversize gaps, incomplete native
bracketing and invalid prospective limits. These are component tests, not a
native Slurm or full-scale overhead experiment.

This closes a raw-evidence gap and provides interval-by-interval process-policy
evaluation. It is not yet wired into final execution admission. Remaining work
includes freezing reviewed bounds, binding the stream review to execution
evidence, checking executable/configuration drift and pressure over the run,
validating the revised native handoff and full-scale overhead, and obtaining a
quiet host. Periodic samples can miss short-lived processes and cannot prove
continuous exclusivity or absence of device contention. No production run was
launched, unrelated process stopped, or DGX accessed.
