# Actual Native-Count Paired Intervals

## Prespecified Successor

`bootstrap_composed_native_qfo_swiss.py` is a prospective executable consumer,
not an executed analysis. It complements the strict retained-interval binder:
when final native counts genuinely differ, old intervals cannot be attached.
This workflow recomputes planned contrasts from the actually admitted native
count audits, without substituting retained values for missing predictions.
Neither passing tests nor preparing this code supplies a final benchmark score.

The existing composed snapshot/audit binding is validated first, including
the separate ordinary/recovered/allocated/composed source routes, raw-file hashes,
reference universe, native admissions and failed-timing distinctions. The new
analysis independently validates confusion counts, family statistics, reference
members/truth totals and aggregates using the unchanged original validation
kernel. Retained counts provide the fixed reference universe for validation
only; unobserved cell values are discarded before draws or contrast calculation.
If no contrast is estimable, refuse to generate unneeded draws.

## Statistical Scope

- Resampling units: all18SwissTrees reference families, not dependent pairs.
- Replicates:100000shared multinomial draws across all observed native cells.
- Seed:20260922 with NumPy PCG64, unchanged from the planned analysis.
- Statistic: mean family PPV and mean family TPR, followed by their harmonic
  mean for F1; recomputed within each replicate. Not mean family F1.
- Contrasts: the unchanged14factorial definitions, including candidate-by-
  reconciliation interactions. Missing required cells retain null metrics.
- Intervals: linear percentile quantiles at0.025/0.975, with planned42endpoint
  Bonferroni quantiles at0.05/84 and1-0.05/84. Never shrink adjustment because
  fewer contrasts are available.
- Report effect sizes, per-family differences and positive/tie/negative counts
  using the original1e-10tie tolerance. Interaction signs are not tool rankings.

These are **new draws from native counts**, not historical interval reuse or
independent confirmation. Report `new_bootstrap_draws=100000`,
`retained_intervals_reused=false` and `unobserved_cells_imputed=false`.
The18development-exposed families retain exchangeability, selection and
approximate percentile-coverage limitations. This supplies no interval for
other QfO challenges, TreeFam-A families or the secondary six-metric mean.

## Execution After Actual Admission

Use the same actual composed score snapshot and actual four source-route audit
path/digest pairs described in
`COMPOSED_NATIVE_SWISS_UNCERTAINTY_PROTOCOL_20261008.md`.
Retain inference24036's successful terminal review, native-pair conversion,
complete six-endpoint scoring and independent admission before using its raw
counts. Failed24038scoring remains missing; no retry is authorized here.

```bash
python -B -m benchmark_tools.bootstrap_composed_native_qfo_swiss \
  --snapshot "$SNAPSHOT" --snapshot-sha256 "$SNAPSHOT_SHA" \
  --retained-counts benchmark_tools/results/qfo_corrected_factorial_complete_20260919/swiss_counts.json \
  --bootstrap benchmark_tools/results/qfo_corrected_factorial_complete_20260919/swiss_bootstrap.json \
  --counts-audit "$ORDINARY_AUDIT" "$ORDINARY_SHA" \
  --counts-audit "$RECOVERED_AUDIT" "$RECOVERED_SHA" \
  --counts-audit "$ALLOCATED_AUDIT" "$ALLOCATED_SHA" \
  --counts-audit "$FINAL_AUDIT" "$FINAL_AUDIT_SHA" \
  --output "$FRESH_NATIVE_INTERVALS"
```

Use the retained scientific Python3.10 environment, disabled user-site/import/
dynamic-loader injection, and one BLAS/OpenMP thread. Existing output paths
are refused. A scheduled execution needs owned submission/release provenance
and its actual resource accounting. Bootstrap analysis is separate from
inference or scoring timing. Shared-host contention is not a blocking gate.

## Validation

Synthetic tests exercise the real validation, aggregate, contrast and sampling
kernels. With all8cells observed, all points, interval metrics and family
differences reproduce the original complete100000draw analysis exactly.
With only native10/12observed, only`C_at_P1_R1`is estimated. A nonuniform
native-count change verifies that macro F1 differs from mean family F1 and
retained F1, with nondegenerate new intervals. Invalid counts, altered reference
members/truth, stale statistics, aggregates, missing families and unestimable
contrasts are rejected. Raw-count handoff tests validate actual changed raw
counts through the separately tested composed admission interface; the
scheduler/admission context is explicitly synthetic.
