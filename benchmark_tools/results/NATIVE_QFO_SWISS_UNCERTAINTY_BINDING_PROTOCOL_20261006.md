# Native SwissTrees Retained-Interval Binding

`benchmark_tools.bind_native_qfo_swiss_uncertainty` links independently admitted
native accuracy and verified family-count audits to the already completed
corrected-release SwissTrees bootstrap. It does not draw another sample,
recount previously audited raw files or create an independent confirmation.

Require the current combined scientific snapshot and its direct replay;
ordinary and recovered count audits have distinct schemas and current sources.
Each audited cell must match its current snapshot's original admission,
native identity and endpoint. Check all recorded source/input/output bindings,
the exact reference and count-conversion convention, and the original bootstrap
count manifest. A recovered count audit must retain null resources and false
timing admission and eligibility.

Require identical full family records, including raw confusion counts,
represented genes and prior-adjusted statistics, plus identical aggregates,
before attaching a retained interval. Aggregate agreement alone is insufficient.
Genuine differing counts remain represented as unmatched; do not overwrite them
or force agreement. Contrasts needing missing or differing cells retain null
metrics and explicit reasons. No historical P1C0R0 admission is invented.

Retain the original 100,000 shared family draws, seed 20260922, percentile
method and adjustment over all 42 planned endpoints. Check the original 14
contrast definitions and recompute their point estimates from retained counts.
Reuse the original intervals and family wins/ties/losses only for fully matched
contrasts; available native comparisons do not define a smaller multiplicity
family. This is SwissTrees only, not uncertainty for VGNC, TreeFam-A, GO, EC,
FAS or the secondary mean.

The previous VGNC rare-error and shared-clade coverage failures remain failures.
This helper does not substitute a new sampling law for them. The original
18-family bootstrap remains development-exposed and conditional on family
exchangeability, with approximate percentile-coverage and selection limits.
It does not establish universal superiority, independent validation, identical
pair decisions or isolated timing. R also changes group-clique to resolved-pair
semantics; P-off retains the initial HMM search.

Validation: 416 joined tests pass in 7.57s, including 36 new binding cases.
Cover normal/recovered audits, actual differing counts, metadata agreement,
duplicate cells, changed provenance/identity/resources, frozen multiplicity,
contrast arithmetic, unavailable cells, interval reuse and no overwrite.
The synthetic test bootstrap is computed once per suite and is not native
benchmark evidence. Prepared Python 3.10 CLI imports pass.

Prospective command after real admission, export and recovered family audit:

```bash
benchmarks/work/native_factorial_review_py310_20261004/bin/python -B \
  -m benchmark_tools.bind_native_qfo_swiss_uncertainty \
  --snapshot ACTUAL_COMBINED_SNAPSHOT/report.json \
  --snapshot-sha256 ACTUAL_REPORT_SHA256 \
  --counts-audit benchmark_tools/results/native_qfo_swiss_family_counts_20261005_v1.json \
    4a49225b50de31f153c1b701f866db7acdfdbefa3d6398339ef5bf8711b2b561 \
  --counts-audit ACTUAL_RECOVERED_COUNT_AUDIT.json ACTUAL_COUNT_AUDIT_SHA256 \
  --retained-counts benchmark_tools/results/qfo_corrected_factorial_complete_20260919/swiss_counts.json \
  --bootstrap benchmark_tools/results/qfo_corrected_factorial_complete_20260919/swiss_bootstrap.json \
  --output FRESH_NATIVE_SWISS_BINDING.json
```

At preparation, scoring22448 remains running and admission22449 pending.
No actual recovered score, recovered count match or native interval binding
is admitted by these fixtures or this command template. Preserve the frozen
scientific protocol and original jobs; do not retry for shared-host contention.
