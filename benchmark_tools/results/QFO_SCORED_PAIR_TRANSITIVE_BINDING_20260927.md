# GO/EC Raw Tables Bound Through Execution Reports

The earlier [direct-record audit](qfo_scored_pair_raw_binding_20260927.json)
looked only for raw files in admissions' direct `checked_records` and
`metric_files` lists. Its zero-match result remains correct for that scope,
but is not the end of the provenance chain.

For six methods, admissions pin an `execution_report`; SonicParanoid and
ProteinOrtho pin their execution results through `checked_records`. Each of
these eight reports contains an `outputs` inventory with the original GO/EC
raw-table length and SHA256. The strengthened runner now requires exactly
one matching raw-output pin and validates execution report bytes before using
it, then rechecks all retained source/input records after analysis.

The [new panel](qfo_scored_pair_panel_20260927_v2.json) successfully verifies
all 16 raw tables against that historical chain. All 56 comparison objects
are exactly equal to the earlier panel. Fifty-one record entries, including
the eight execution reports, were verified. No raw file, native endpoint,
benchmark default or previous receipt was changed. No native scoring or
inference was repeated.

This resolves the specific missing raw-byte binding for this panel. It does
not imply complete transitive provenance of every input, immutable external
attestation, correct annotation semantics or independent confidence intervals.

## Reproduction

```bash
python -m benchmark_tools.run_qfo_scored_pair_panel --repo . \
  --output NEW_BOUND_PANEL.json
python -m pytest -q tests/unit/test_run_qfo_scored_pair_panel.py \
  tests/unit/test_compare_qfo_scored_pairs.py tests/unit/test_audit_qfo_go_ec.py
```

All 37 focused tests pass. New coverage rejects altered pair identities even
when counts and means are unchanged, and exercises the checked-record route
to execution reports. These additions postdate the latest full unit suite;
no new full-suite result is implied.
