# Three-Run Native Archive Audit

Added `audit_dual_native_diagnostics.py` to compose the provenance, raw CPU
replay, native-output validation and canonical-output comparison checks.
All three assigned jobs must have detailed terminal scheduler records before
the driver reads any native archive. A failed task remains in the result.

The complete archived recipe is checked against its pinned inventory. For
each successful task, every regular file in its run directory is hashed
before inspection and checked again afterward; symlinks and inventory changes
are rejected. Checked native products and input evidence are also retained.
This inventory proof does not imply semantic validation of every intermediate
tree or alignment: native semantics are checked by `validate_scaling_outputs`.

Canonical outputs are compared with prior periodic tasks 1, 3 and 8, as
prespecified. The prior complete pressure-panel audit is SHA-256-pinned to
`d656ae37dcb64d617a132c82391745e33bfd5f078f6c1df7186156f0022ea006`;
the selected prior task evidence is rehashed during comparison. The actual
retained audit digest and all three validated comparison entries were checked
while preparing this driver.

The report retains native wall time, GNU-time accounting, cgroup memory,
original and narrow CPU flags, complete screening/pressure results, output
counts, output identities, file inventories and temporal-order issues.
Output mismatch is a failure with observations preserved. CPU flags are
reported independently of file/output validation, and no result admits
scientific timings, environmental isolation or publication readiness.

## Validation

All 96 focused tests passed in 1.44 seconds across the audit, provenance,
raw replay, launcher, collector and dual-reader modules. The orchestration
fixtures mock the expensive native-output parsing and raw replay; their
validators have separate focused tests. Actual complete archive validation
remains pending until all three diagnostic jobs terminate.

Invocation after collection:

```sh
python -m benchmark_tools.audit_dual_native_diagnostics \
  --archive benchmarks/work/dual_native_archive_21912 \
  --results benchmark_tools/results \
  --scheduler-directory benchmarks/work/dual_native_scheduler_21912 \
  --prior-audit benchmark_tools/results/dgx_pressure_overhead_audit_21889_20260919.json.gz \
  --output benchmarks/work/dual_native_audit_21912_v1.json
```

The output must not already exist. Preserve failures and raw evidence; do
not retry selectively or substitute later runs for this prospective set.
This three-run diagnostic is not the repeated overhead panel and does not
fulfill the 27 matched-resource scaling runs required for publication.
