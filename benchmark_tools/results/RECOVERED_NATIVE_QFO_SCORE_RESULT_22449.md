# Recovered Native QfO Accuracy Result

Original native 22437 produced recoverable P0/C0/R1 predictions, but its
allocation/timing failed. Conversion 22447, fresh scoring 22448 and independent
admission 22449 now complete. Scoring elapsed 32:09 on 8 CPU/64 GiB; admission
elapsed 0:24 on 2 CPU/32 GiB, both COMPLETED 0:0. These are postprocessing
allocations, not repaired or valid native inference timing.

The [published admission](recovered_native_qfo_assessment_admission_22449.json)
matches the original 781,476-byte report, SHA
`4c3a17a76eff5043c8b40d0c9f1e8ead6c7ed32d4e248dbd485988c33537ae3b`.
All 15 tasks are fresh/noncached and 48 native metric records checked;
1,702 records bind admission. Original inference is not restarted and
FAILED 1:0/null/ineligible timing persists. No default or endpoint changes.

## Admitted Scores

Both ablations retain initial HMM search; downstream profile refinement (P)
and candidate expansion (C) are off. R changes group-clique predictions to
the native inferred-pair strategy, not just a clustering parameter. These are
development-exposed ablations, not selected tool defaults.

| Metric | P0/C0/R0 | P0/C0/R1 |
| --- | ---: | ---: |
| VGNC F1 | 0.666834 | 0.898185 |
| SwissTrees F1 | 0.689184 | 0.789574 |
| TreeFam-A F1 | 0.605404 | 0.602508 |
| GO similarity | 0.472119 | 0.490260 |
| EC similarity | 0.932114 | 0.967702 |
| FAS | 0.774945 | 0.785007 |
| Secondary six-metric mean | 0.690100 | 0.755539 |
| Submitted pairs | 9,009,082 | 5,113,820 |
| Inputs with any relation | 55.5258% | 55.1107% |

[Actual seven-row export](native_qfo_scientific_scores_20261006_v1/scores.md),
[full-precision TSV](native_qfo_scientific_scores_20261006_v1/scores.tsv) and
[bound snapshot](native_qfo_scientific_scores_20261006_v1/report.json) preserve
P0/C0/R0. Five other fresh native score rows remain unavailable, not zero, live
status or substituted historical scores. All-input relation coverage uses
984,137 inference proteins and is not accuracy. GO/EC/FAS are not F1; the
six-metric mean is a project-defined secondary summary.

P0/C0/R1 precision/recall are 0.99953910/0.81549260 for VGNC,
0.94915197/0.67593220 for SwissTrees and 0.95668249/0.43971892 for TreeFam-A.
The small TreeFam-A point decrease is retained; gains are not uniform.

## SwissTrees Uncertainty

The [family audit](recovered_native_qfo_swiss_counts_22449_20261006.json)
checks all 18 families and 10,765 scored reference relations. Full records,
including counts, members and prior-adjusted statistics, match the retained
corrected cell. No normal native raw count audit is repeated.

The [guarded binding](native_qfo_swiss_uncertainty_binding_22449_20261006.json)
matches two cells and only R-at-P0/C0. It reuses the original 100,000 shared
family draws, seed 20260922, preserving all 42 adjusted endpoints. The other
13 native contrasts remain unavailable. No new bootstrap is drawn.

| R-on minus R-off | Difference (pp) | Adjusted interval (pp) | Positive/tie/negative families |
| --- | ---: | --- | --- |
| F1 | +10.0390 | [-4.6418, 24.4199] | 14/1/3 |
| Precision | +30.5212 | [14.3252, 46.5473] | 16/1/1 |
| Recall | -6.5334 | [-18.0413, -0.0436] | 0/10/8 |

This supports a conditional precision-recall trade-off, not a clear adjusted
F1 gain, independent confirmation or universal advantage. Family exchangeability,
development exposure and approximate percentile-coverage limits remain.
Count-derived F1 uses full precision and does not replace rounded endpoints.
No interval for other QfO challenges or the secondary mean is admitted;
VGNC rare-error/shared-clade uncertainty failures remain unresolved.

## Independent Readback

[Score readback](recovered_native_qfo_score_readback_22449_20261006.json) checks
original/published admission, actual scheduler outcomes, 48 raw metric records,
15 task rows, six aggregations, rational-decimal arithmetic and all 112 TSV
fields. This is bounded readback, not another full transitive raw admission
or independent biological validation.

[SwissTrees readback](recovered_native_qfo_swiss_readback_22449_20261006.json)
independently accumulates new raw counts/members, checks 10,765 pairs and 18
families, exact record matches, rational aggregate/contrast arithmetic, family
signs and retained intervals. Maximum family arithmetic discrepancy is
2.22e-16. Reference-tree truth reconstruction is not repeated; exact pair-label
anchor checks remain in the primary audit.

The FAS sample contains 252,451 pairs and 257,077 proteins (4.9366% of eligible
pairs); 126,521 proteins occur in multiple sampled pairs, maximum degree 17.
Native pair-IID SEM 0.0003342373 is not a paired-method or family CI. Sampling
is unseeded; the 9000 cap concerns new missing-pair computations, not total
sampled pairs. Missing annotations/scores and population limits remain.

Validation before use: 416 joined tests, including 36 new binding cases, plus
12 guide contracts pass. Post-result integrated validation passes 436 tests
in 9.53s, including 8 new portable artifact checks. Those check retained
report/table/document consistency, not raw admission or biological validity.
All 920 frozen helpers and 19 plan evidence files
remain unchanged. Native 22444 remains live; 22445 is its original dependent
review. No next native identity before terminal review/fresh gates. Full
scientific, uncertainty, provenance, manuscript and release requirements
remain; this result does not establish publication readiness.

Timing measurements were collected on a shared Threadripper while other
analyses were running. Competition for CPU, memory bandwidth and I/O may have
affected elapsed times, with an unknown and potentially tool-dependent impact.
These are observed shared-host timings, not estimates of isolated performance.
