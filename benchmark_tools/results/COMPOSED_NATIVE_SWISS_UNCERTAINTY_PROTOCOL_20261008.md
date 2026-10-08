# Composed-Native SwissTrees Uncertainty

## Scope

This is a prospective consumer, not an executed count audit or a new accuracy
admission. Native inference job 24036 must finish successfully, pass its prepared
terminal review, and complete native-pair conversion, all six QfO assessments
and independent admission before a final-cell snapshot or raw audit can be used.
Native11 scoring job 24038 failed with OUT_OF_MEMORY; its five completed endpoint
tasks supply no admitted score, family audit or contrast. Do not retry it here.

The new count audit accepts only integer index12/p1_c1_r1 and the explicit
`composed_native_qfo_scientific_reporting_v1` snapshot. It replays the composed
reporter with the saved final admission, verifies bound outputs and sources,
then reads only the selected checked, inventoried SwissTrees raw file. It uses
the unchanged original reference/count/native-endpoint kernels and the retained
18-family reference universe. F1 must reproduce the admitted SwissTrees endpoint
within the existing absolute tolerance of 5e-8.

The successor binder accepts ordinary, recovered, allocated and composed audits
as distinct source/schema routes. Prior cell audits are reused, not recounted.
Their admission, endpoint, reference, raw-file and source bindings must agree
with the new snapshot. Recovered scientific accuracy does not restore failed
timing eligibility or resource measurements.

## Interval Rules

Reuse the unchanged original interval projector only when every family record
AND aggregate agrees exactly with the retained counts. Aggregate equality alone,
even with unchanged family confusion counts but different represented members,
does not establish a match. Missing or differing cells produce null contrast
metrics, never cached substitutions. The full 14 contrasts and 42 metric
endpoints remain in the multiplicity adjustment, even when only some are
estimable. The retained seed20260922, 100000 shared draws, alpha0.05 and linear
quantiles cannot change. These consumers perform zero new bootstrap draws.

If final native counts differ, retain that difference. A genuinely new paired
analysis must be explicitly implemented and validated rather than presenting
the old intervals as applicable. No uncertainty is supplied for GO, EC, FAS,
VGNC, TreeFam-A or the project-defined secondary six-metric mean by this workflow.
SwissTrees remains development-exposed, with only18families, approximate
percentile coverage and family-exchangeability limitations. This is not
independent confirmation or a publication-readiness claim.

## Commands After Successful Admission

Use a fresh composed score export and its actual digest. Set `SNAPSHOT` and
`SNAPSHOT_SHA` to that export, not synthetic test data. Set `FINAL_AUDIT` to a
fresh result path. Preserve all original files and failed results.

```bash
python -B -m benchmark_tools.audit_composed_native_qfo_swiss_counts \
  --snapshot "$SNAPSHOT" --snapshot-sha256 "$SNAPSHOT_SHA" \
  --retained-counts benchmark_tools/results/qfo_corrected_factorial_complete_20260919/swiss_counts.json \
  --index 12 --output "$FINAL_AUDIT"
```

The retained counts digest is
`c9d8bf02ef6f287c56d073fa61f165c982caf834a94fc6e44e6589f03e46eba6`.
The retained bootstrap digest is
`f777dead1294b3810c0479877aade7fb4411c0544696f671226ff198432e9211`.
The binder pins both and verifies the bootstrap's original count binding.

Supply actual path/digest pairs for the prior ordinary, recovered and allocated
audits, and the new final audit. Their retained current paths are respectively:

- `benchmark_tools/results/native_qfo_swiss_family_counts_20261005_v1.json`
- `benchmark_tools/results/recovered_native_qfo_swiss_counts_22449_20261006.json`
- `benchmark_tools/results/native10_allocated_swiss_counts_20261007_v1.json`

```bash
python -B -m benchmark_tools.bind_composed_native_qfo_swiss_uncertainty \
  --snapshot "$SNAPSHOT" --snapshot-sha256 "$SNAPSHOT_SHA" \
  --retained-counts benchmark_tools/results/qfo_corrected_factorial_complete_20260919/swiss_counts.json \
  --bootstrap benchmark_tools/results/qfo_corrected_factorial_complete_20260919/swiss_bootstrap.json \
  --counts-audit "$ORDINARY_AUDIT" "$ORDINARY_SHA" \
  --counts-audit "$RECOVERED_AUDIT" "$RECOVERED_SHA" \
  --counts-audit "$ALLOCATED_AUDIT" "$ALLOCATED_SHA" \
  --counts-audit "$FINAL_AUDIT" "$FINAL_AUDIT_SHA" \
  --output "$FRESH_BINDING"
```

Use the retained scientific Python3.10 environment with user-site imports and
Python/dynamic-loader injection disabled; constrain BLAS/OpenMP to one thread.
Any scheduled execution must retain owned submission/release provenance and
validate the actual resource envelope without requiring an uncontended host.
Count reconstruction, interval binding, inference and QfO scoring are distinct
resource scopes. No inference or scoring job is launched by these commands.

## Test Interpretation

Tests explicitly stub the admission/reporting handoff, not the raw count or
interval kernels. They use real18-family reconstruction and100000draw paired
bootstrap fixtures. A supported mixed allocated/composed contrast is
`C_at_P1_R1`; its metrics and family differences must reproduce retained values
exactly. Other missing contrasts remain null. Production reporter/native
admission contracts are tested separately. Passing synthetic tests is not a
final native benchmark result.
