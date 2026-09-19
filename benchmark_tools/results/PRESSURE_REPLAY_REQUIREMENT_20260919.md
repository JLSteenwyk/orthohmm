# Pressure Replay Requirement

Prepared while array21889 and its controller recorder21890 remain live.
No DGX files, logs or partial native outcomes were accessed. This changes
only the local post-run evaluator, not the frozen deployed collectors,
native workload, pressure diagnostic interpretation or overhead budgets.

`replay_frontier_overhead_measurement.replay` now accepts the keyword-only
`expected_native_pressure` argument:

- `True` requires native pressure evidence in every raw observation point.
- `False` rejects pressure evidence in a legacy-only expected collector.
- `None` preserves the existing automatic replay behavior for old callers.

The completed pressure-panel auditor must pass `True`. Without this explicit
expectation, deleting pressure from every observation and replacing the
retained screening with a consistent legacy screening could be accepted as
a legacy measurement. Partial pressure coverage was already rejected by the
collectors' evaluators; the new gate also rejects complete absence when
pressure was required. Integer/string lookalikes are not boolean settings.

After the coverage gate, existing replay still requires raw/report agreement,
recomputes CPU and pressure diagnostics, validates identities, scopes,
counter totals and memory evidence, and requires exact screening equality.
No numerical tolerances, exclusions or environmental-admission rules changed.
Boundary measurements continue to lack periodic interval coverage.

## Validation

```sh
/home/bizon/anaconda3/bin/python -m pytest -q \
  tests/unit/test_replay_pressure_overhead.py \
  tests/unit/test_replay_frontier_overhead_measurement.py \
  tests/unit/test_frontier_pressure_integration.py \
  tests/unit/test_audit_frontier_overhead.py
```

Result:89passed in7.48seconds. The24new cases use synthetic pressure-enabled
boundary/periodic fixtures and cover successful recomputation, missing/partial
pressure, unexpected pressure, invalid expectation types, altered reported
deltas, raw totals, cgroup identities and scopes. Historical replay and
complete-panel audit tests still pass. These are evaluator tests, not results
from21889 or evidence that it meets the engineering budget.

## Remaining Panel Integration

The old complete-panel auditor is still bound to21838's plan, authorization,
recipe paths and scheduler command. Do not invoke it unchanged for21889.
Before auditing the new panel, bind its exact submitted identities and call
the strict pressure replay gate for every candidate successful task. Preserve
the all18terminal gate before native inspection, detailed scheduler evidence,
native output equivalence, duration/budget checks and every failed task.
Exact same-runtime replay remains required; prior cross-runtime floating-point
differences must not be silently relaxed. No scientific timing or publication
readiness is admitted by this prerequisite.
