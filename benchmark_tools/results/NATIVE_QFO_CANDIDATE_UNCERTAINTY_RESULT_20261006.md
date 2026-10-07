# Native Candidate-Expansion Uncertainty

## Actual New Evidence

The third admitted fresh native QfO cell is P0C1R0, not a live-job or cached
score substitution. Its six-endpoint results remain in the unchanged
[current three-cell snapshot](native_qfo_scientific_scores_20261006_v2/scores.md).
The earlier scientific admissions and recovered failed-timing row remain
unchanged. No native inference, conversion, assessment or FAS sampling rerun.

Use the original exporter to supply only the new candidate admission to the
original ordinary count auditor. This deliberately isolated
[audit input](native_qfo_candidate_count_input_20261006_v1/report.json) is
not a replacement comparison table: other admitted cells are omitted to avoid
recounting their completed audits. Its report SHA256 is
`5bd9eaeea790a54b472275b14d33bddfa267f7f4923d57b0fff98baa54c215d5`.
The source and original scientific kernels remain unchanged.

The [new count audit](native_qfo_candidate_swiss_counts_20261006_v1.json)
reads the actual inventoried SwissTrees raw output of assessment 23894.
All 18 complete family records, represented members and aggregate exactly
match the retained P0C1R0 records. It preserves native decimal F1
0.6857094982951092; count-derived full-precision F1 is 0.6857094995663944.
This tiny serialization difference does not define a new endpoint.

The unchanged [guarded binder](native_qfo_candidate_swiss_uncertainty_20261006_v1.json)
combines the new audit with the existing P0C0R0 and recovered P0C0R1 audits.
Exactly two contrasts now have matched native records: C_at_P0_R0 and
R_at_P0_C0. All other contrasts remain unavailable with null metrics.
Reconciliation metrics remain exactly the retained values. Reuse the frozen
100,000 draws and seed 20260922, retaining all 42 endpoints in adjustment;
no new draws, uncertainty for other QfO endpoints or selection-adjusted
independent validation is created.

## Independent Readback

The [new stdlib readback](native_qfo_candidate_swiss_readback_20261006_v1.json)
uses separate CSV/gzip parsing and Fraction arithmetic, not the primary
family-count or bootstrap functions. It checks 36 native family records and
32,295 raw rows: native baseline, native candidate and retained candidate.
All 10,765 candidate scored pair labels, truth and represented members match
the retained evidence. Reconstruct macro precision/recall first and their
harmonic mean second; do not average family F1. Independently reproduce the
candidate effect and family wins/ties/losses. Intervals are checked against
the frozen bootstrap output, not independently recalculated or certified.
This is not a transitive native admission or whole-partition equivalence.

An initial real readback exits 1 with `Changed bootstrap provenance` before
creating a result file: the original bootstrap count reference points to
the byte-identical work-directory copy, while the retained copy is under
results. Update only the new readback to check both exact records, byte length,
digest and parsed content; retain the original data, auditor and binder.
The final real readback exits 0. Tests cover this distinct-path equality and
changed-original, changed-supplied, different-byte and parsed-content failures.
The initial 96-pass test receipt predates this recovery; the expanded suite
passes 101 tests in 4.62s with zero failures/errors/skips.
The final [joined suite](native_qfo_candidate_uncertainty_joined_tests_20261006_v1.xml)
passes 256 tests in 7.26s, including exact manuscript score/interval tables
and current claim-scope contracts. Saved readback is 53,823 bytes, SHA256
`07941ba9c55ba7dc14c5ba53df37d28cb3b2f60f599b9d177b10c9afa99c1329`.
Direct readback confirms the final source still matches and the entire prior
R_at_P0_C0 contrast is unchanged, with exactly two matched and 12 unavailable
contrasts. These tests/readbacks do not resolve family exchangeability.

## Conditional Result

| Metric | Difference (pp) | Nominal 95% Paired CI | 42-Endpoint Adjusted CI | Family Wins/Ties/Losses |
| --- | ---: | --- | --- | --- |
| SwissTrees F1 | -0.347446 | [-3.332692, 1.955055] | [-5.816832, 3.309051] | 5/8/5 |
| SwissTrees precision | -3.523453 | [-7.426534, -0.291442] | [-10.229612, 0.373979] | 4/8/6 |
| SwissTrees recall | +4.375165 | [2.100638, 6.872093] | [0.940111, 8.688718] | 9/9/0 |

Candidate expansion raises recall, with a lower precision point estimate;
the adjusted F1 interval includes zero. This does not prove equivalent F1
or general superiority. All 18 families remain development-exposed, with
exchangeability and approximate percentile coverage limits. P-off retains
initial HMM search; both R-off configurations submit group-clique pairs.
The mixed other-endpoint results, unseeded FAS limitations and secondary
mean interpretation remain in the admitted-result document. No default tuning.

## Reproduce Readback

Run the checked existing audits/binder inputs or the final independent
readback with a fresh output path; existing paths refuse overwrite.
The retained outputs already exist and need not be regenerated on continuation.

```bash
env -u PYTHONPATH -u PYTHONHOME -u PYTHONUSERBASE \
  -u LD_PRELOAD -u LD_LIBRARY_PATH -u LD_AUDIT \
  PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  benchmarks/work/native_factorial_review_py310_20261004/bin/python -B \
  -m benchmark_tools.readback_native_qfo_candidate_swiss \
  --audit benchmark_tools/results/native_qfo_candidate_swiss_counts_20261006_v1.json \
  --audit-sha256 7ed2ea1a7a94d55b28ecbadd7bd2be8e99df64b9fab3c266ebd987e8dbcfc479 \
  --binding benchmark_tools/results/native_qfo_candidate_swiss_uncertainty_20261006_v1.json \
  --binding-sha256 76e6f4d33a13a12b7f0a396f6b795b63b654666a038485fb8af2538259f553c2 \
  --output benchmark_tools/results/native_qfo_candidate_swiss_readback_20261006_v1.json
```

The source/test/report milestone updates the manuscript and claims checklist,
not the existing two-cell figures, supplement or final release payloads.
Four fresh score identities remain unavailable; identity 9's pre-native
failure is retained and the same identity-10 job 23902 remains live at this
checkpoint. Observe that handle without resubmission; 11/12 stay sequential
behind reviewed history. Full scientific and publication requirements remain.

Timing measurements were collected on a shared Threadripper while other
analyses were running. Competition for CPU, memory bandwidth and I/O may have
affected elapsed times, with an unknown and potentially tool-dependent impact.
These are observed shared-host timings, not estimates of isolated performance.
