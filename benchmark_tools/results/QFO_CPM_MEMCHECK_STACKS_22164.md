# Memcheck 22164 Stack Inventory

Extracted the retained 166,213,125-byte XML using a streaming standard-library
parser, checking its SHA256 before and after extraction:
`ebf7a23098fcf3f0135083370f783ad010116760ce4c99dab5f3d9c1ee4b53db`.
The [machine-readable inventory](qfo_cpm_memcheck_stacks_22164.json) preserves
all non-leak stacks, auxiliary allocation text, error occurrence counts and
leak top-frame counts. Ten unit tests cover extraction, separate record/event
counts, malformed/incomplete evidence, identity changes and fatal-signal retention.

## Observations

- Four InvalidRead and five UninitCondition records, each with one reported
  occurrence. All nine primary stacks have `__wcscat_avx2` at the top and
  include `getpath.c` and `Py_InitializeFromConfig` frames.
- These observed stacks locate reports in libc wide-string operations called
  during the retained Python executable's initialization. They do not point
  to the later `read_partition` garbage-collection crash site recorded for
  admission 22155.
- The XML also contains 8,128 Leak_PossiblyLost and 62 Leak_DefinitelyLost
  records. The allocation wrapper at the top of a leak stack does not identify
  the responsible owner. Leak stacks remain in the original XML; this summary
  neither dismisses them nor establishes their relation to the crash.
- The XML records RUNNING then FINISHED and no fatal_signal element. The
  instrumented child still returned 97; this report does not change that failure.

## Reproduce

```bash
python -m benchmark_tools.summarize_cpm_memcheck \
  --xml benchmarks/results/qfo_cpm_refinement_memcheck_diagnostic_v1/memcheck.xml \
  --sha256 ebf7a23098fcf3f0135083370f783ad010116760ce4c99dab5f3d9c1ee4b53db \
  --output /tmp/qfo-cpm-memcheck-stacks-new.json
python -m pytest -q tests/unit/test_summarize_cpm_memcheck.py
```

## Interpretation and Next Check

The nine non-leak reports are not evidence that the OrthoHMM refinement code
caused the earlier SIGSEGV. Conversely, startup-location evidence does not
establish false positives, memory safety, or a repair. The original scientific
admission remains failed and the high-CPM accuracy endpoint remains missing.

The next bounded diagnostic should compare one minimal, no-site-import Python
startup (`-S -c pass`) under the same interpreter and Memcheck instrumentation
to these recorded stacks, retaining even a negative reproduction. This would
test whether any reports reproduce without loading OrthoHMM or its scientific
extensions; it would not test refinement correctness or justify bypassing
admission. No such control was run in this step. No scientific workload was
repeated. Post-failure input/runtime/partition revalidation remains outstanding.
