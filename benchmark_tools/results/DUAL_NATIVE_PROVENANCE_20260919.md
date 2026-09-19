# Dual Native Provenance Checks

Added `verify_dual_native_provenance.py` without changing the frozen launcher,
collector or running scheduler recorder. This complements raw measurement
replay; neither checker alone admits scientific timings.

The context loader pins the plan and recipe and reconstructs the plan from
the pinned prior pressure panel and prospective protocol. The verifier binds
each task to its assigned non-array job (21912, 21913, 21914), successful
terminal state, exclusive 20-CPU/96-GiB DGX allocation, zero restarts and exact
batch-script path. Duplicate scheduler fields and array records are rejected.

Preparation must match the planned native command, GNU-time wrapper, input
checksums, enumeration order and, for OrthoFinder, copied-input paths/order.
Before/after runtime records must match the two frozen runtime manifests and
the execution recipe. OrthoFinder additionally requires post-run copied-input
verification. The wrapper source, embedded measurement and launched collector
must match. Receipt and runtime comparisons distinguish JSON booleans from
integers; finite positive preparation/check durations are required.

Validation:

```sh
pytest -q tests/unit/test_verify_dual_native_provenance.py tests/unit/test_replay_dual_native_measurement.py tests/unit/test_dual_native_diagnostic.py tests/unit/test_measure_native_dual_bracket_step.py tests/unit/test_probe_dual_cpu_brackets.py
```

All 84 tests passed in 1.06 seconds. The first actual retained terminal record,
`benchmarks/work/dual_native_scheduler_21912/scheduler_21912.txt`, passed the
scheduler identity check for diagnostic index zero. No remote native records
were inspected while the diagnostic quiet window remained active.

## Remaining Checks

After all three jobs terminate, collect immutable archives and apply these
checks to actual preparation/verification/task/measurement records. Validate
the complete archived recipe and input bytes, replay raw observations, check
native artifact inventories and compare canonical native outputs with prior
periodic tasks 1, 3 and 8. Report all failures, interval flags, pressure totals,
duration and memory without selective reruns or retrospective policy changes.

Retained before/after runtime records cannot rule out temporary changes during
execution. Exact metadata matching is not proof of native output correctness,
environmental isolation, collector overhead bounds or publication readiness.
