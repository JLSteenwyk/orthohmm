# Five Native Cells Support A Bounded Profile Contrast

The [five-result native snapshot](native_factorial_progress_20261005_v5/report.json)
adds terminal-reviewed/scoredP1C0R1. The unchanged binding helper checks
all70 RefOG records per cell, including metadata, against the pinned
eight-cell factorial. All350 records and F1/P/R estimates match exactly.
The [new binding](native_factorial_uncertainty_binding_20261005_v3.json)
therefore supports five conditional contrasts, retaining the original12
planned contrasts and36-endpoint multiplicity adjustment. Seven native
contrasts remain unavailable; their missing-cell lists reflect the new input.
No native score is imputed and older snapshots/bindings remain unchanged.

New contrast: downstream profile expansion on versus off, with candidate
expansion off and reconciliation on (P1C0R1 minusP0C0R1).

| Metric | Difference (pp) | Nominal 95% Paired CI | Adjusted Percentile CI |
| --- | ---: | --- | --- |
| F1 | 0.6063680296846456 | [-0.7818660834460829, 2.829213190422542] | [-1.2114359037526952, 4.520970311951387] |
| Precision | -0.03625742104762253 | [-0.8130543786057117, 0.66817855793141] | [-1.3160359491892883, 1.144899575220177] |
| Recall | 0.9115348176624067 | [-0.7756114915307013, 3.892888286648632] | [-1.1841218051888265, 6.317027679625983] |

Family F1 wins/ties/losses are4/61/5. Every nominal and adjusted interval
for this contrast spans zero. These are reused deterministic intervals from
the original20,000 paired RefOG draws/seed20260918/alpha0.05, not new draws
or independent confirmation. Original weighted F1 is recomputed within
each draw rather than averaging family F1. The prior four supported
contrasts retain their exact numerical results.

Both cells retain initial sensitive HMM search. This measures a conditional
downstream profile effect, not the total HMM contribution or superiority
over OrthoFinder. Root-HOG co-membership is not resolved ortholog-pair truth.
Development exposure, possible nonexchangeable families and approximate
percentile/simultaneous coverage remain. Failed22427 retains its scientific
recovery, not successful timing. No timing estimate or correction is made.

The [independent readback](native_factorial_uncertainty_readback_20261005.json)
checks all350 records and independently recomputes weighted F1/P/R within
1e-10 percentage points. It validates exact retained intervals, family
counts, failed-wrapper status, missing native cells and920 unchanged helpers.
This is direct evidence binding, not another transitive raw review.

The extended manuscript and claim-to-evidence table now include this neutral
profile result and the distinction between five native-backed contrasts
and the full retained eight-cell analysis. No figure, PDF, archive, default
or frozen scientific/helper source is replaced. Full-goal requirements remain.

## Validation

All24 joined tests pass, zero failures/errors/skips,0.65s: four new document
bindings,16 existing uncertainty projection controls and four existing
candidate-manuscript controls. The [JUnit](native_uncertainty_manuscript_tests_20261005.xml)
is retained. Tests verify the exact rounded effect/interval, zero-spanning
interpretation,70-family counts, missing native scope, source links and
independent readback. Existing controls reject changed family metadata or
counts despite equal aggregates. After editing, all920 native helper pins
and the queued reviewer source remain unchanged. Fresh Slurm confirms22433
RUNNING10:07 with actualPID4028693 and22434PENDING; preserve those handles.

## Reproduce

Use a fresh output path and the original retained local evidence:

```bash
benchmarks/work/release_alert_refresh_20261001/venv/bin/python -B \
  -m benchmark_tools.bind_native_orthobench_uncertainty \
  --snapshot benchmark_tools/results/native_factorial_progress_20261005_v5/report.json \
  --snapshot-sha256 46c6d5f79cbac94b90d717e2cfc7f9b4d352f899f557291759e7e43d70f33705 \
  --factorial benchmark_tools/results/orthobench_factorial_results_20260916.json \
  --output /tmp/native_orthobench_uncertainty_fresh_v3.json
```

This does not certify portable raw data, distribution rights or publication
readiness. Preserve actual live22433 and dependent22434; do not rerun or
duplicate them to produce this statistical summary.
