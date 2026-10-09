# Native QfO Terminal-Failure Addendum

### Terminal Native QfO Outcomes

The later terminal reviews do not add any admitted scores: four of seven
fresh cells retain admitted accuracy and three retain missing scores.
The preceding four-cell endpoint table, figure, intervals and localization
remain unchanged. This update is failure reporting, not a complete factorial,
independent confirmation or superiority over OrthoFinder.

| Cell | Accuracy status | Six-metric mean |
| --- | --- | ---: |
| p0_c0_r0 | supplied_native_admission | 0.69009982 |
| p0_c0_r1 | supplied_recovered_scientific_admission | 0.75553924 |
| p0_c1_r0 | supplied_native_admission | 0.68471435 |
| p0_c1_r1 | no_supplied_native_admission | Unavailable |
| p1_c0_r1 | supplied_allocated_native_admission | 0.75463772 |
| p1_c1_r0 | scoring_OUT_OF_MEMORY_no_admission | Unavailable |
| p1_c1_r1 | native_SIGSEGV_no_admission | Unavailable |

Native11 (P1/C1/R0) completed inference and group-clique conversion,
retaining 11,755,521 pairs and relation coverage 585,180/984,137.
Its assessment ended OUT_OF_MEMORY at 8 CPU slots/32 GiB; FAS exited137.
Five other endpoint tasks completed but supply no admitted score or mean.
Coverage is not accuracy; its exact scoring peak memory remains unknown.

Native12 (P1/C1/R1; job24036) ended FAILED1:0. The observer native-step
wrapper completed0:0, but the actual native command exited-11 (SIGSEGV).
The final native receipt remained running and final metrics/orthogroup
outputs were absent. Retained downstream phylogeny files are partial
evidence, not successful inference or admitted pair predictions. Conversion
and scoring were not run. No stack trace was captured; the root cause and
crash location are unknown. The last buffered progress message does not
identify the crashing operation. No automatic inference retry was performed.

Reviewed failed-command observations: 64526.103655 wall seconds,
1739733.457142 native-task CPU seconds and 20492017664 bytes
native-step lifetime peak memory, with the documented wrapper/launcher scopes.
These describe a failed attempt, not successful inference or scoring cost.

Shared Threadripper measurements while other analyses were running. CPU, memory-bandwidth and I/O contention may affect elapsed times by an unknown, potentially tool-dependent amount. These are observed shared-host timings, not estimates of isolated performance.

Missing outcomes remain null, not zero. The three supported SwissTrees
contrasts and all 42 planned adjusted endpoints retain their original scope;
no new bootstrap draws, endpoints, functional similarity claims or default
selection follow from this update. Prepared success-only five-cell reporting
was not executed. This revision remains not submission-ready.

Evidence:
- [snapshot](../native_qfo_scientific_scores_20261007_v3/report.json) (SHA256 `7aab1cbd31fb650df42e6e80a14ee0167b35a810f08202768c89c514755211a2`)
- [native11](../native11_qfo_scoring_failure_addendum_20261008_v1/report.json) (SHA256 `4f1a5c536af8ce992d46866374b3ab44ee3158077316d69ba2e9f5fb34cee203`)
- [review](../../../benchmarks/work/native12_composed_terminal_review_20261008_v1/review.json) (SHA256 `a55a6da5ffb8dfbb362883937c915b161796a9ddd7275f4ac48a1c4d049e5c4a`)
