# Preserve Failed Snapshot Evidence

Commit `bf4e102` changes failure reporting, not counter admission. The
frontier sampler now retains its inventories and successfully read counter
prefix in a `FrontierSnapshotError` when topology changes or a counter
disappears during observation. It still raises an error. Periodic and
boundary collectors save that evidence and the preceding hierarchy
observation to `failed_point.json` before rethrowing. Existing cleanup still
releases the owned worker; no successful measurement report is emitted.

The artifact is explicitly marked invalid and not admitted for scientific
timing. A missing after-inventory is null, not reconstructed. No retry chooses
a more favorable snapshot, missing counters are not zero-filled, and the
original comparison and CPU screening thresholds are unchanged. This does
not repair task 6 of array 21838: its missing within-snapshot inventories
remain missing. The new code has not been deployed into an active frozen
benchmark recipe or used to admit timing results.

## Verification

Focused tests initially passed 95 cases, including within-snapshot topology
change, disappearance during a counter read, saved evidence before rethrow,
refusal to overwrite a failure artifact, worker cleanup, and historical
numerical replay. The full suite exposed two historical-source identity
tests that expected current collector files to remain byte-identical forever.

Preserved the three old source files from committed revision `9e89534` as
test fixtures. Historical recipe tests continue checking their original
byte counts and SHA-256 hashes, not newly invented hashes. Current-source
behavior tests and historical numerical replay remain separate. The initial
full-run result (two failures) remains in
`benchmarks/work/frontier_failure_diagnostics_unit_20260918.xml`.

After those test corrections, the final full suite passed **5,711 tests with
9 skips in 122.89 seconds**:

```bash
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /home/bizon/anaconda3/bin/python -m pytest -q tests/unit \
  --junitxml=benchmarks/work/frontier_failure_diagnostics_unit_final_20260918.xml
```

Final JUnit SHA-256:
`4f0d6b63a1daf87558519f9f254a2abb98b0f0f2feaade69bfddef2e41cdd2f1`.
No live timing experiment or native CLI rerun was needed for this
diagnostic-only change. A prospective resource-measurement design that
handles transient services without losing their activity is still required.
