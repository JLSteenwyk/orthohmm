# GDB Diagnostic 22159: Crash Not Reproduced

The single prespecified refinement-only GDB execution completed in 44 seconds
with exit 0:0, using one CPU and 64 GiB. The child exited normally and wrote its
completion JSON. No crash stack was obtained; the post-exit register request
reported that no registers were available. The retained shared-library inventory
includes libraries without debugging information and a linux-vdso symbol warning.

Rechecked all 986 recorded input/source references, child log and completion
report. Independently repeated the output/metadata check against the successful
native refinement: 984137 unique genes in 390845 groups, identical 23875927-byte
partition with SHA256
`f4c6f1973bc9636828baf1e6d9be3f416a18fc8ada5302fb502c082495fc1811`.
Runtime verification before and after the diagnostic agrees.

[Submission](qfo_cpm_backtrace_submission_22159.json) and
[checked result](qfo_cpm_backtrace_result_22159.json) retain exact executor,
debugger, scheduler and artifact identities. The executor is the clean detached
`cpm_backtrace_diagnostic_v1_20260923` checkout at
`4a414d59ac59f4a8a74550758bab47a1ed4b8e84`; its 21 focused tests passed before
submission, including installed-debugger failure/signal smoke tests.

This negative reproduction does not identify the crash cause or prove the
runtime is safe. Instrumentation can change memory layout and timing. Do not
repeat until a desired crash or successful result appears. Failed admission
22155 remains failed, candidate 22156 remains blocked, and the high-CPM accuracy
endpoint remains missing. No optimization, scoring, parameter/default change or
DGX action occurred. A smaller independent reproducer or concrete native-memory
evidence is still needed before claiming a fix or revising scientific admission.
