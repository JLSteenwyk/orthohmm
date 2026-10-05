# Full-Native Factorial Progress

Timing measurements were collected on a shared Threadripper while other analyses were running. Competition for CPU, memory bandwidth and I/O may have affected elapsed times, with an unknown and potentially tool-dependent impact. These are observed shared-host timings, not estimates of isolated performance.

P toggles downstream profiles, C candidate expansion, R reconciliation. Initial HMM search remains on.
OrthoBench F1 is reference-group co-membership, not native resolved-pair accuracy.

| Index | Dataset | Cell | Outcome | F1 (%) | Wall (s) | CPU (s) | Peak (GiB) |
| ---: | --- | --- | --- | ---: | ---: | ---: | ---: |
| 0 | orthobench | p0_c0_r0 | failed_wrapper_science_recovered | 69.7634 | 2289.7647 | 67271.6775 | 4.6140 |
| 1 | orthobench | p0_c0_r1 | native_success | 72.7050 | 3254.9227 | 93796.6952 | 4.6299 |
| 2 | orthobench | p0_c1_r0 | native_success | 66.6395 | 2306.6354 | 67315.2179 | 4.6278 |
| 3 | orthobench | p0_c1_r1 | native_success | 73.4023 | 3522.0292 | 102817.7650 | 4.7245 |
| 4 | orthobench | p1_c0_r1 | native_success | 73.3114 | 3460.6570 | 97852.1129 | 10.5875 |
| 5 | orthobench | p1_c1_r0 | no_supplied_terminal_review | Unavailable | Unavailable | Unavailable | Unavailable |
| 6 | qfo_corrected | p0_c0_r0 | no_supplied_terminal_review | Unavailable | Unavailable | Unavailable | Unavailable |
| 7 | qfo_corrected | p0_c0_r1 | no_supplied_terminal_review | Unavailable | Unavailable | Unavailable | Unavailable |
| 8 | qfo_corrected | p0_c1_r0 | no_supplied_terminal_review | Unavailable | Unavailable | Unavailable | Unavailable |
| 9 | qfo_corrected | p0_c1_r1 | no_supplied_terminal_review | Unavailable | Unavailable | Unavailable | Unavailable |
| 10 | qfo_corrected | p1_c0_r1 | no_supplied_terminal_review | Unavailable | Unavailable | Unavailable | Unavailable |
| 11 | qfo_corrected | p1_c1_r0 | no_supplied_terminal_review | Unavailable | Unavailable | Unavailable | Unavailable |
| 12 | qfo_corrected | p1_c1_r1 | no_supplied_terminal_review | Unavailable | Unavailable | Unavailable | Unavailable |

CPU includes the wrapper; peak includes the native-step launcher, not pure algorithm RSS.
Preparation, conversion and scoring are outside the native interval.

- Only explicitly supplied terminal reviews are summarized; blanks are not live-job status or zeros.
- Recovered failed-wrapper resources are failed-attempt observations, not clean-success timings.
- Missing accuracy is not inherited from cached predictions; QfO metrics require separate native admission.
- This checks direct report identities, not transitive raw artifacts or new independent accuracy.
- No averages, causal component overhead, isolated ranking, retries or job release.
- The three reused configurations and historical cached stage costs remain separately reported.
