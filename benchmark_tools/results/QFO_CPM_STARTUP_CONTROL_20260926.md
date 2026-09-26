# Minimal Python Startup Control

One local diagnostic ran the retained interpreter with `-B -S -c pass` under
the same Memcheck flags as diagnostic 22164. No site initialization, OrthoHMM
code or scientific-package imports were requested. The runner reused the
allocator and thread overrides, but used the checkpoint-recovery directory
for both cwd and PYTHONPATH instead of the original replay directory. This
recorded difference means the control is not an exact environment match.
It checked the retained
Python/libc/Valgrind component identities before and after execution, imposed
a 120-second timeout, and did not retry. This was not a scientific workload
or a matched-resource measurement.

The [retained report](qfo_cpm_startup_control_20260926.json) records exit 97,
normal XML RUNNING/FINISHED states, zero fatal-signal records and 3.51 seconds
of descriptive wall time. The XML digest is
`649d8b283c6e9912a1a14e4ac6148f72eca07bb694eda80c57cf87dc8a515893`
(14,755,353 bytes, retained locally rather than committed).

| XML error records | Refinement diagnostic 22164 | Minimal startup |
| --- | ---: | ---: |
| InvalidRead | 4 | 8 |
| UninitCondition | 5 | 3 |
| Leak_PossiblyLost | 8,128 | 2,017 |
| Leak_DefinitelyLost | 62 | 0 |

All 11 startup-control non-leak reports have `__wcscat_avx2` top frames and
include `getpath.c` and `Py_InitializeFromConfig`, as do all nine retained
refinement non-leak reports. This is a reproduction of error kinds and stack
locations, **not identical stacks, counts, allocations or a reproduction of
the later SIGSEGV**. The control shows that loading OrthoHMM/scientific
extensions is not necessary to obtain these startup reports on this runtime.
It does not establish false positives or rule out independent refinement bugs.
Differences in leak counts are not causal effects or proof of ownership.

Fourteen focused runner/extractor tests pass, including one-attempt execution,
nonzero outcome retention, timeout evidence, source identity, malformed XML,
and refusal to reuse the output directory.

## Reproduction Command

Executed once from the repository root:

```bash
python -m benchmark_tools.run_cpm_startup_control \
  --status benchmarks/results/qfo_cpm_refinement_memcheck_diagnostic_v1/status.json \
  --output benchmarks/work/qfo_cpm_startup_control_20260926
```

The runner intentionally refuses that existing directory. The raw diagnostic
and interpreter installation are external requirements; this is not a portable
runtime-validation bundle. No suppression, scientific parameter change or
revised admission was introduced. High-CPM scores remain missing, and the
original failure remains failed. The earlier crash still requires a targeted
reproducer or native crash evidence; repeated successful runs cannot replace it.
