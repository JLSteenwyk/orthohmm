# Native Reviewer Signal-11 Diagnostic

## Retained Failure And Scope

On goal resumption, accounting confirms native22444 COMPLETED0:0, elapsed
11:42:59 (native step11:18:09), reviewer22445 FAILED0:11, elapsed31:25.
Conversion22450, assessment22451 and admission22452 each FAILED1:0 at their
predecessor-accounting gates. Preserve these outcomes and every original file.
No successful review.json or failure.json exists in the original review
directory; runtime, resource replay, resources and environment intermediates
exist, but do not constitute a completed review or scientific admission.

Kernel journal records python PID818632 segmentation fault at15:20:23.
addr2line on the retained Python binary maps the instruction pointer0x53ea83
to _Py_INCREF. This is a crash location, not an explanation of its cause.
No matching Python core was found in /var/crash; coredumpctl is unavailable.
The original reviewer stderr is empty. Do not infer an OOM, a NumPy fault,
input corruption or algorithm failure from these observations alone.

## One Bounded Diagnostic

Run the unchanged original validate_native_factorial_outputs production CLI
alone under the original Python3.10 review venv, with Python fault reporting
and unbuffered output enabled. Use a new Slurm job on bizon/gpu with2 CPUs,
32GiB RAM,6hours, no requeue. Check fresh safe capacity and unchanged request,
plan, source and helper bindings before release. Do not modify services,
unrelated jobs, frozen scientific sources, inputs or existing result namespaces.

Pinned request SHA256:
35f8ef9b1d7abc1f574b8e7e3c55bc91c68c9f8f74247a9a816754a41a93fb6a.
Pinned validator SHA256:
3357637503f35238f654edea5c4c12bd293f8d20908c141421a2d40278af7c8d.
Pinned original reviewer SHA256:
63e7d7fdda52afa7a36eecd89492a2684e8310260febbfb5adf6d23f88fd26de.
Pinned plan SHA256:
6c87babcbb5581830e0b9e7b9bf9aaba30a85bde1c4ab465e561017e67e9c89b.

Write only new diagnostic outputs under
benchmarks/work/native_review_sigsegv_diagnostic_20261006_v1. The semantic
CLI performs its original terminal, runtime/source/input/partition/checkpoint
and candidate gates. Record /usr/bin/time, stdout, stderr and scheduler outcome.
This diagnostic does not rerun native inference, resource replay, conversion,
scoring or admission. It is not an automatic retry of the failed reviewer.

If it fails, retain all evidence and use its stack or error to choose the next
action. If it succeeds, conclude only that standalone semantic validation
succeeds in this invocation; do not claim the original crash is explained,
fixed, harmless or a passing full review. Consider a separately documented
postprocessing repair that preserves all original terminal-review gates.
Do not release native index9 or scoring using this diagnostic alone.

## Reporting

All elapsed times here are shared-host diagnostic observations, not inference
timings or isolated performance. Competing CPU, memory-bandwidth and I/O
demands may affect them by an unknown, potentially tool-dependent amount.
The full publication goal remains active and completion unproven.
