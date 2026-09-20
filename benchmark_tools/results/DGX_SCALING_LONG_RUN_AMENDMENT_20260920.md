# Long-Run Collector Amendment Before Execution

The unexecuted replacement plan v1 selected the historical root-context
wrapper, which only accepts 900 seconds. Its underlying worker also rejects
timeouts above 900 seconds, and historical replay only accepts 60/900 seconds.
The previous adapter tests intercepted measurement and did not exercise this
incompatible boundary. No replacement native run has been launched.

Select the separate `measure_scaling_root_context.measure` entry point in v2.
It accepts exactly 20 CPUs, 96 GiB, an 85,800-second native timeout and a
one-second requested counter interval. It reuses the existing point readers,
CPU/pressure evaluators, owned-process-group timeout cleanup and memory reads,
but provides a long-run worker/observer lifecycle. Six-digit point filenames
preserve chronological lexical ordering above 9,999 observations. If the
initial observation fails, the worker receives an abort gate, not permission
to launch an unobserved native command. Native failure/timeout outcomes remain.

No native method, input, task identity, output path, allocation, timer boundary
or scientific setting changes. The original v1 specification remains retained.
The amendment does not authorize execution, service changes or timing admission.

Historical replay must not be relaxed in place. A new matching long-run replay
is required before launch, as are composed lifecycle tests and a fresh source
recipe. The prior four-proteome overhead results validate the earlier collector
configuration, not the new long-run lifecycle or its worst-case memory footprint.
Points are retained in memory and evaluated after execution; this cost must
remain explicit. Larger-input and long-run behavior need validation.

The existing point readers observe host counters, not an exhaustive process,
GPU or device-I/O inventory. Passing `monitor_host=True` through the generic
measurement API does not create such an inventory. The environment-policy
and whole-run workload evidence requirements therefore remain unresolved.
Do not treat a parameter name as proof of monitoring coverage.
