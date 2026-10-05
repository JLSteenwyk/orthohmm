# Six Native Cells Support Two Bounded Profile Contrasts

The [six-result native snapshot](native_factorial_progress_20261005_v6/report.json)
adds terminal-reviewed/scoredP1C1R0. The unchanged binding helper checks
all70 RefOG records per cell, including metadata, against the pinned
eight-cell factorial. All420 records and F1/P/R estimates match exactly.
The [new binding](native_factorial_uncertainty_binding_20261005_v4.json)
supports six conditional contrasts, retaining the original12 planned
contrasts and36-endpoint multiplicity adjustment. Six native contrasts
remain unavailable. No score is imputed; older snapshots and bindings remain.

New contrast: downstream profile expansion on versus off, with candidate
expansion on and reconciliation off (P1C1R0 minusP0C1R0).

| Metric | Difference (pp) | Nominal 95% Paired CI | Adjusted Percentile CI |
| --- | ---: | --- | --- |
| F1 | 0.5679664267983782 | [-0.6952861922180759, 2.650746205783999] | [-1.1132088102248157, 4.290615410522566] |
| Precision | 0.2118681102550255 | [-0.6972810715641037, 1.4185920927778046] | [-1.1512970442046537, 2.252320179955322] |
| Recall | 0.977755327473318 | [-0.7453284646915769, 4.0025729294105785] | [-1.1577634680346018, 6.4922540494016285] |

Family F1 wins/ties/losses are6/59/5. Every nominal and adjusted interval
for this new contrast spans zero. These intervals are reused from the
original20,000 paired RefOG draws/seed20260918/alpha0.05, not new draws or
independent confirmation. Original weighted F1 is recomputed within each
draw, not averaged across family F1 values. The prior five supported
contrasts retain their exact numerical results, including the C0R1 profile
effect+0.606368pp with adjusted interval[-1.211436,4.520970].

Both cells retain initial sensitive HMM search. This measures a conditional
downstream profile effect, not the total HMM contribution or superiority
over OrthoFinder. Final group co-membership is not resolved ortholog-pair
truth. Development exposure, possible nonexchangeable families and approximate
percentile/simultaneous coverage remain. Failed22427's scientific recovery
does not become successful timing. No timing estimate or correction is made.
The60 changed non-reference genes inP1C1R0 do not change its RefOG records;
this binding does not infer the cause of those partition differences.

The [independent readback](native_factorial_uncertainty_readback_20261005_v2.json)
checks420 records and independently recomputes weighted F1/P/R within1e-10
percentage points. It validates exact retained intervals, family counts,
the five unchanged previous contrasts, exact missing-cell lists, failed-
wrapper status and920 unchanged helpers. This is direct evidence binding,
not another transitive raw review or new bootstrap.

The extended manuscript and claim-to-evidence table now include the new
neutral profile result and distinguish six native-backed contrasts from
the full retained eight-cell analysis. No figure, PDF, archive, default or
frozen scientific/helper source is replaced. Full-goal requirements remain.
Preserve the original liveQfO22435; no inference/scoring/launch occurs in
this statistical update.

## Validation

All25 joined tests pass, zero failures/errors/skips,0.66s: five manuscript
evidence tests,16 uncertainty projection controls and four candidate-
manuscript controls. The [new JUnit](native_uncertainty_manuscript_tests_20261005_v2.xml)
is retained separately. Checks cover both rounded profile effects, their
zero-spanning interpretation, fixed factor settings, family counts, missing
native scope, source links and independent readback. Existing projection
controls reject changed family metadata/counts despite equal aggregate scores.
After edits all920 native helper pins remain unchanged, as does the7,210-byte
binding helper SHAff6c6f3ac394f930044c7108352d508bf42f90bd322e61c4b004545e4edb3fe4.
Fresh scheduler confirms22435RUNNING8:42 with actualPID4059977, matching
creation time and affinity0..31. Preserve that original live attempt.

## Reproduce

Use a fresh output path and the original retained local evidence:

```bash
benchmarks/work/release_alert_refresh_20261001/venv/bin/python -B \
  -m benchmark_tools.bind_native_orthobench_uncertainty \
  --snapshot benchmark_tools/results/native_factorial_progress_20261005_v6/report.json \
  --snapshot-sha256 4024277cdfcb39d109b37360fb2efc09be93af5cdbbb83cfffdccf88a1315a1f \
  --factorial benchmark_tools/results/orthobench_factorial_results_20260916.json \
  --output /tmp/native_orthobench_uncertainty_fresh_v4.json
```

This does not certify portable raw data, distribution rights or publication
readiness. Retained older bindings and JUnit artifacts are not overwritten.
