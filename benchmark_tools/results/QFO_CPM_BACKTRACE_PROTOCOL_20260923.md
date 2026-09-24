# One Native Backtrace Diagnostic

The allocator-debug refinement 22158 failed with SIGSEGV after writing the
partition, without an allocator violation report or completion JSON. No Python
core is present in the retained crash directory. The next diagnostic is exactly
one refinement-only child under GDB, using the same saved graph/profile inputs,
frozen scientific runner, runtime, one-thread settings, PYTHONMALLOC=debug and
PYTHONFAULTHANDLER=1 as 22158. No optimizer, threshold, seed or garbage-collection
setting is changed. It is not a retry of scientific admission.

Allocate one CPU, 64 GiB, at most one hour on bizon, with no requeue. Use a fresh
`qfo_cpm_refinement_backtrace_diagnostic_v1` directory. Check input/runtime/source
identities before launching; record GDB binary identity and version. Disable
GDB init files, automatic script loading and debuginfod network downloads.
Preserve ordinary address randomization. Run once and request all-thread native
backtraces, registers and loaded shared-library information at termination.
Capture combined debugger/child output and the debugger return status.

A signal stop, debugger failure, missing child completion report or changed
input must remain a failed diagnostic, never a scientific success. If the child
completes, retain the negative reproduction outcome and independently compare
its output and metadata; success does not repair failed admission 22155. Do not
retry until a crash or successful partition is obtained. Do not release candidate
job 22156. The complete high-CPM accuracy endpoint remains missing.

GDB can affect timing and memory layout. A backtrace localizes a crash, not
necessarily the corrupting operation or root cause. No timing comparison or
memory-safety claim follows. Any code fix or revised scientific admission needs
separate review and validation. DGX remains deferred.
