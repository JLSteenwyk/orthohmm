# Boundary Arm For Lineage Overhead Measurement

`measure_lineage_boundary_step.py` prepares the boundary-only arm needed for
a future complete overhead panel. It imports the exact same point reader,
native worker and memory reader as the periodic lineage collector. Its
measurement loop retains command release, timer, one-second completion polling
and cleanup, but reads counters only before release and after completion.
No periodic sample is taken while the native command runs.

It writes `lineage_boundary_report.json` with `lineage_boundary_screening_v1`.
The two points use the same `native_lineage_v1` format, including mandatory
pressure and both host endpoints. Evaluation requires exactly two points
enclosing the complete native command, retains the same whole-command CPU
screen, and reports hierarchy, lineage and pressure changes. Interval
screening is explicitly unavailable and interval flags are null, not an
empty list interpreted as a quiet run. This arm measures neither a monitor-free
baseline nor cold-cache execution.

`replay_lineage_boundary_measurement.py` uses the same command/job/native-time,
raw/report equality, memory, evidence-hash and admission checks as the periodic
lineage replay. Its inventory requires exactly two point files: zero and a
positive canonical final polling index. Noncontiguous indices are expected
because completion polls do not produce samples. It recomputes the boundary
screen without inventing periodic coverage; failures and changed files reject
replay. Scheduler, runtime and native output provenance remain separate audits.

## Verification

All 90 targeted tests pass, including 27 new boundary/replay tests:

```sh
python -m pytest -q tests/unit/test_measure_lineage_boundary_step.py \
  tests/unit/test_replay_lineage_boundary_measurement.py \
  tests/unit/test_measure_native_lineage_step.py \
  tests/unit/test_replay_lineage_native_measurement.py \
  tests/unit/test_measure_frontier_boundary_step.py
```

Tests check exact worker/reader sharing and lifecycle correspondence, only
two reads despite multiple completion polls, success/failure/timeout cleanup,
whole-screen equality, required pressure and identity consistency, missing or
extra samples, invalid poll indices, changed commands and memory, altered raw
records, false admission and changes during replay. The fixtures use synthetic
lineage additions and mock workers; there is no live boundary-arm result yet.

## Deployment Status

This code was developed locally while native diagnostic jobs 21995-21997 ran.
It was not copied into their frozen deployment and did not change their
collector, protocol or endpoints. No overhead task has been submitted using
this boundary arm. Wait for all current diagnostics to finish, collect and
audit their complete evidence, then freeze the complete future overhead panel
with both lineage arms. Existing failed panels and flags remain unchanged.
Scientific timing admission and publication readiness are not established.
