# Complete-Command Native Pressure Integration

Added opt-in `native_pressure=True` to the periodic frontier collector and
its boundary-only control. The default remains false. Previously exported
recipes and historical reports are not modified or reclassified.

Each enabled observation now includes native-step CPU, memory and I/O PSI
after the cgroup frontier read and before the final outer host snapshot.
The existing probe checks raw fields, process membership, cgroup identity
and per-read timestamps. Integration additionally verifies matching scope,
clock ticks, boot identity and enclosure within the complete observation.
Periodic evaluation reports adjacent-interval and whole-command deltas;
boundary evaluation reports only the whole-command delta and continues to
declare interval screening unavailable.

Partial pressure coverage, changing native cgroup identity, counter resets,
inconsistent raw fields and invalid timing brackets fail evaluation. A
pressure read failure retains the preceding hierarchy/frontier observation
and error in `failed_point.json`, then propagates to the existing worker
cleanup. It does not fabricate the missing pressure fields.

Original CPU screening and admission flags remain unchanged. Pressure is
diagnostic, not a new exclusion threshold, subtraction-based estimate of
foreign work, or wall-time correction. Native wrappers/descendants contribute
to the measured scope. Added reads can increase collector overhead.

## Verification

```bash
python3 -m pytest -q \
  tests/unit/test_frontier_pressure_integration.py \
  tests/unit/test_measure_native_frontier_step.py \
  tests/unit/test_probe_native_pressure.py \
  tests/unit/test_measure_frontier_boundary_step.py \
  tests/unit/test_verify_frontier_overhead_provenance.py \
  tests/unit/test_replay_frontier_overhead_measurement.py \
  tests/unit/test_run_dgx_frontier_overhead.py
```

167 tests passed in8.14s. The17 new integration cases cover bracketed reads,
scope/identity/order/counter faults, mixed coverage, retained failures,
whole-command versus interval reporting, and both collector lifecycles at
native exit0,7,124. Source-contract tests preserve worker lifecycle and the
boundary/periodic distinction. Historical evidence tests still use retained
source fixtures; their hashes and scoring rules were not changed.

During development, one lifecycle-source comparison needed updating for the
new optional reader. Three added boundary lifecycle tests initially failed
because the shared mock setup expected an unused `read_hierarchy` attribute;
the test adapter now follows the existing boundary-test pattern. These were
test integration failures, not successful live runs.

## Remaining Work

This change has not yet been deployed as a complete-command DGX run.
Export the new collectors with their transitive source dependencies,
including `probe_native_pressure.py` and `audit_dgx_pressure.py`, to a fresh
directory. Verify a short complete command before freezing another overhead
panel and scientific inclusion rules. The old eighteen-run panel remains
incomplete/invalid for its intended overhead conclusion, and the original
27 scaling runs remain unadmitted. Neither unit tests nor the prior short
CPU-pressure injection controls establish controlled inference timings.
