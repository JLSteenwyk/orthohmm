# Native Lineage Archive Audit Preparation

## Frozen Bindings

`verify_lineage_native_provenance.py` binds only corrected jobs 21995,
21996 and 21997 to their frozen three-task plan and recipe. It reconstructs
the plan from the pinned baseline and protocol, verifies receipt/task/argv,
input checksums and order, OrthoFinder copied-input checks, before/after
runtime records, wrapper identity, worker/collector launch and scheduler
resources. It requires exclusive spark-7ff0, 20 CPUs/96 GiB, no restart or
requeue, successful completion and exact non-array job identity.

The scheduler `Command` is the actual local submission-script path, not the
remote deployed recipe path. Tests explicitly distinguish these. The failed
launch jobs 21990-21992 cannot substitute for the corrected jobs. Plan and
recipe hashes remain those in `LINEAGE_NATIVE_SUBMISSION_21995.md`; no deployed
source, native command, protocol or threshold has been changed by this work.

## Archive Checks

`audit_lineage_native_diagnostics.py` requires all three detailed terminal
scheduler records before it reads the native archive. It validates the
archived recipe inventory against the pinned manifest, hashes each native
run directory before and after processing, binds provenance, invokes the
lineage raw-measurement replay, validates native products and GNU-time
evidence, and computes canonical output fingerprints.

The prior same-method comparison is the existing pressure-panel v2 audit:
`dgx_pressure_overhead_audit_21889_20260919.json.gz`, SHA-256
`d656ae37dcb64d617a132c82391745e33bfd5f078f6c1df7186156f0022ea006`.
This checksum was reverified locally. Comparison uses its original periodic
task indices 1, 3 and 8, retaining their evidence checks. Canonical equivalence
does not prove biological accuracy or identical internal work.

All three assigned tasks remain in the result, including scheduler failures,
missing reports, invalid measurements and output mismatches. Original and
narrow CPU flags remain explicit. Missing tasks produce null aggregate
conclusions rather than a successful subset. Temporal order and boot-domain
discrepancies are reported separately. Scientific timing admission,
environmental validity and publication readiness remain false.

After all three jobs and controller capture are terminal, collect a complete
archive preserving paths relative to the DGX project root, including the
deployed recipe, native output roots and plan-referenced input files. Keep the
controller scheduler directory separate. Then use a fresh output path:

```sh
python -B -m benchmark_tools.audit_lineage_native_diagnostics \
  --archive benchmarks/work/lineage_native_archive_21995 \
  --results benchmark_tools/results \
  --scheduler-directory benchmarks/work/lineage_native_scheduler_21995 \
  --prior-audit benchmark_tools/results/dgx_pressure_overhead_audit_21889_20260919.json.gz \
  --output benchmarks/work/lineage_native_audit_21995.json
```

This command has not been run against the current native outputs: the panel
is still active, and no DGX output inspection has occurred during it.

## Tests And Limits

All 140 targeted tests pass, including 48 new provenance/archive tests:

```sh
python -m pytest -q tests/unit/test_verify_lineage_native_provenance.py \
  tests/unit/test_audit_lineage_native_diagnostics.py \
  tests/unit/test_replay_lineage_native_measurement.py \
  tests/unit/test_measure_native_lineage_step.py \
  tests/unit/test_verify_dual_native_provenance.py \
  tests/unit/test_audit_dual_native_diagnostics.py
```

New tests cover exact plan/recipe bytes, all three task bindings, changed
input/order/runtime/worker/scheduler fields, failed-attempt substitution,
premature native inspection, missing reports, altered archives, replay
failure, output mismatch, retained flags and changed clock domains.
Archive-assembly tests mock native product validation and raw replay; their
unit results are not evidence that the real jobs pass those checks. Existing
native-product validators remain separately tested and will run on collection.
The retained runtime records do not exclude temporary runtime changes, and
the archive audit does not establish overhead or absence of interference.
