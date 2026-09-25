# One Refinement-Only Memcheck Diagnostic

The ordinary and allocator-debug refinement admissions failed with SIGSEGV;
the one-shot GDB diagnostic 22159 did not reproduce the crash. No retained
Python core is available. Run exactly one saved-input refinement child under
the installed Valgrind Memcheck to seek native-memory evidence. This is not
another scientific admission attempt and cannot release candidate 22156.

Use the frozen scientific runner and identical saved graph/profile inputs,
one CPU, 64 GiB RAM, one-hour limit, no requeue, and a fresh
`qfo_cpm_refinement_memcheck_diagnostic_v1` directory. Do not rerun the optimizer,
change clustering parameters, score output or retry until a preferred outcome.
Verify source/input/runtime identities before and after execution, retain
Valgrind binary/version identity, raw XML and combined child output.

Use `PYTHONMALLOC=malloc` and `PYTHONFAULTHANDLER=1`, with one-thread library
settings. Disabling Python's small-object allocator is an instrumentation
change, not a scientific-method change or a repair. The
[CPython Valgrind guide](https://github.com/python/cpython/blob/main/Misc/README.valgrind)
explains the allocator concern; the
[Valgrind manual](https://valgrind.org/docs/manual/manual-core.html) documents
XML reporting and nonzero error-exit status. Record origins and up to 40 stack
frames. Ignore user configuration, do not trace exec children, and use no added
suppression file. Leak detection is outside this diagnostic's scope; memory
access and undefined-value reports remain enabled. Use error exit code 97.

A nonzero child/tool exit, missing or malformed XML, incomplete XML run state,
reported memory error, changed inputs, missing child completion metadata or
changed output bytes must remain a failed diagnostic. A clean result requires
normal completion and the same independent partition/metadata checks used in
22159, but still does not prove memory safety or repair admission 22155.
Instrumentation changes timing and memory layout. Reports may originate in
Python, libraries or OrthoHMM; classify actual stacks before attributing cause.
Do not infer comparative runtime, algorithmic superiority or publication readiness.
