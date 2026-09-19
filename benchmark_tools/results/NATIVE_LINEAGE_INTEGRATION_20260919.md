# Native Lineage Collector Integration

## Scope

`measure_native_lineage_step.py` provides a separate complete-command
engineering collector. It uses the existing worker, Slurm step lifecycle,
native timer and final cgroup memory reader. The worker lifecycle is copied
from the dual-bracket collector with only report-name and description changes;
a source comparison test checks that correspondence. No frozen existing
collector or CPU threshold file is changed.

The new format identifies points as `native_lineage_v1`, screening as
`native_lineage_screening_v1` and interval results as
`native_lineage_interval_v1`. It writes `lineage_report.json`; the new
`replay_lineage_native_measurement.py` produces
`lineage_native_measurement_replayed`. The validator rejects frontier fields
and requires lineage plus pressure evidence. No old frontier failure is
converted or silently retried using the new reader.

## Checks And Retained Evidence

The sequence is native/job hierarchy and original host bracket, aggregate
root-to-job lineage, native pressure, then the extended outer host endpoint.
The saved original endpoint supplies the narrow window; both windows use
the same native CPU sample and timestamps. The validator requires native/job
scope agreement, boot agreement, ordered and enclosed reads, nondecreasing
parent CPU, pressure scope agreement and stable lineage/pressure identities
between samples. Membership is rechecked after observation.

Every failed observation can retain its partial native, hierarchy and lineage
evidence in `failed_point.json`; the worker is released through the existing
cleanup path. Native nonzero exits and timeouts remain failures. Both CPU
screens retain their original thresholds and flagged interval indices.
Lineage aggregate deltas and signed differences are additional diagnostics,
not admission rules. There is no subtraction from native wall time.

Replay requires contiguous raw samples, agreement between embedded and raw
native/memory evidence, expected job/argv/resources/cadence/timeout, zero
native exit without timeout, exact wall-duration recomputation, and matching
screening recomputation. It validates memory scope, observation time and
integer current/peak counters, hashes evidence before and after, and rejects
an altered sample inventory. Replay uses the collector's evaluator, so it
is not an independent statistical implementation or an audit of scheduler,
runtime binaries, scientific output identity or input provenance.

## Verification And Limits

All 118 targeted tests pass, including 48 new integration/replay tests:

```sh
python -m pytest -q tests/unit/test_measure_native_lineage_step.py \
  tests/unit/test_replay_lineage_native_measurement.py \
  tests/unit/test_measure_native_dual_bracket_step.py \
  tests/unit/test_replay_dual_native_measurement.py \
  tests/unit/test_probe_cgroup_lineage.py \
  tests/unit/test_run_lineage_lifecycle_control.py
```

These cover successful, failed and timed-out mock workers, cleanup and failure
retention, timing/identity/counter corruption, schema separation, raw/report
disagreement, command/resource/cadence mismatch, changed files during replay,
and memory errors. Integration fixtures derive from retained hierarchy
observations but use synthetic lineage/pressure additions. They are not
native DGX integration results.

The separate completed-service control 21989 is still the only live DGX
evidence for the new aggregate reader. Next freeze and run an end-to-end
native measurement control for this integrated collector, preserve all raw
outputs and replay them, and evaluate service churn during observation.
Then prepare a complete fresh boundary-versus-periodic overhead panel using
the new schema in both arms and frozen native commands. Neither successful
unit tests nor prior lifecycle results authorize comparative scientific timing.
Full-node behavior, overhead, non-CPU isolation and a prospective scientific
inclusion policy remain open; the publication goal is incomplete.
