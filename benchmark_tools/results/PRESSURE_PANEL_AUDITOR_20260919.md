# Pressure Panel Auditor

Prepared locally while21889remains active, without reading any DGX files or
partial native outcomes. The preceding33e738c pressure-coverage gate is now
used by the complete-panel auditor when selecting `--panel pressure_21889`.
The default `frontier_21838` retains the historical panel binding.

## Frozen Identity

The allowlisted panel specifications bind separate exact plan, recipe and
authorization SHA-256 values, array identities, recipe roots, manifest paths
and scheduler commands. Unknown panels are rejected. The new scheduler
command is the absolute local submission path retained in
[the original21889submission](dgx_pressure_overhead_submission_21889.txt),
not an assumed path inside the DGX recipe. Both prospective pressure
protocols must match their entries in the pinned recipe manifest.

Successful tasks still require exact task receipts, native preparation/input
order, before/after runtime verification, wrapper identity, collector worker
launch, exclusive20CPU/96GiB scheduling with no restart, raw measurement
replay, native-output validity and canonical output fingerprints. Every
pressure-panel successful replay explicitly requires native pressure at all
observation points. Whole-command and periodic-interval pressure diagnostics
are retained separately from unchanged CPU screens. Boundary intervals remain
unavailable; pressure is diagnostic only, not a new exclusion threshold.

All18scheduler tasks must be terminal before context/native archive reads.
Every failed, missing or invalid task remains in the audit. Existing paired
statistics, >=60-second minimum native duration, <=5% per-method median and
<=10% every-pair budgets are unchanged. Numerical budget, observation/output
validity and environmental validity remain separate. Neither this auditor
nor a numerical pass authorizes the scientific27-run scaling experiment.

## Post-Terminal Invocation

Only after local controller accounting establishes all18terminal states:
collect the exact archived recipe and outputs, native raw observations and
the recorder's detailed terminal `scheduler_INDEX.txt` files without rewriting
their contents. Retain accounting and validate recorder completeness. A new
complete audit can then be invoked with paths to those retained artifacts:

```sh
python -m benchmark_tools.audit_frontier_overhead \
  --panel pressure_21889 \
  --archive ARCHIVE_ROOT \
  --results benchmark_tools/results \
  --accounting ACCOUNTING_FILE \
  --output FRESH_AUDIT_JSON
```

These are placeholders, not a record of a completed audit. Exact replay in
the collector's runtime remains required; cross-runtime floating-point
differences must be preserved and explained, not hidden with a new tolerance.
No transfer, actual-panel audit or scientific timing admission occurred here.

## Tests

```sh
/home/bizon/anaconda3/bin/python -m pytest -q \
  tests/unit/test_verify_frontier_overhead_provenance.py \
  tests/unit/test_audit_frontier_overhead.py \
  tests/unit/test_replay_pressure_overhead.py \
  tests/unit/test_replay_frontier_overhead_measurement.py \
  tests/unit/test_summarize_frontier_overhead.py \
  tests/unit/test_pressure_overhead_panel.py \
  tests/unit/test_run_dgx_frontier_overhead.py
```

286passed in8.07seconds. Tests bind all18tasks for each allowlisted panel,
reject cross-panel records and changed checksums/protocols, exercise both
recipe roots against missing/extra/changed/symlink files, ensure terminal
gating precedes native/context reads, verify pressure is mandatory in the
new audit path, and retain failed-task/budget behavior. New full-panel and
provenance fixtures are synthetic; these tests are not observations from the
running experiment. Earlier historical fixture tests remain included.
