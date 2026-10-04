# Retained Factorial Stage Costs

P: profile expansion; C: candidate expansion; R: reconciliation. P-off retains initial HMM search.
All values are historical shared-host observations, not isolated efficiency comparisons.
Preparation repeats the same shared arm in two rows; it excludes search/profile/seed construction.
NA is unavailable/not applicable, not zero. RSS is sampled summed process-tree RSS, not lifetime peak.

| Dataset | Cell | Shared preparation (s) | Reconciliation (s) | Mean CPU cores | Sampled tree RSS (GiB) | Full pipeline |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| OrthoBench | p0_c0_r0 | 0.908 | NA | NA | NA | NA |
| OrthoBench | p0_c0_r1 | 0.908 | 1582.806 | 23.125 | 1.467 | NA |
| OrthoBench | p0_c1_r0 | 14.539 | NA | NA | NA | NA |
| OrthoBench | p0_c1_r1 | 14.539 | 1952.122 | 26.745 | 1.573 | NA |
| OrthoBench | p1_c0_r0 | 0.654 | NA | NA | NA | NA |
| OrthoBench | p1_c0_r1 | 0.654 | 1514.539 | 25.508 | 1.476 | NA |
| OrthoBench | p1_c1_r0 | 11.491 | NA | NA | NA | NA |
| OrthoBench | p1_c1_r1 | 11.491 | 1963.790 | 28.675 | 1.577 | NA |
| Corrected QfO | p0_c0_r0 | 3.452 | NA | NA | NA | NA |
| Corrected QfO | p0_c0_r1 | 3.452 | 4445.429 | 30.643 | 5.925 | NA |
| Corrected QfO | p0_c1_r0 | 79.772 | NA | NA | NA | NA |
| Corrected QfO | p0_c1_r1 | 79.772 | 6861.549 | 30.151 | 7.332 | NA |
| Corrected QfO | p1_c0_r0 | 3.036 | NA | NA | NA | NA |
| Corrected QfO | p1_c0_r1 | 3.036 | 4416.658 | 30.466 | 5.921 | NA |
| Corrected QfO | p1_c1_r0 | 77.641 | NA | NA | NA | NA |
| Corrected QfO | p1_c1_r1 | 77.641 | 6822.365 | 29.871 | 7.308 | NA |

## Limits

- Historical shared-host stage observations, not isolated tool speed or causal component overhead.
- Preparation intervals omit cache loading and seed/profile construction. An arm is shared by its R-off/R-on rows; do not sum them as independent preparation runs.
- Eight R-off reconciliation stages are not applicable, not zero-cost full pipelines.
- Peak is sampled sum of process-tree RSS, potentially double-counting shared pages and missing between-sample or short-lived peaks. It is not the newer scaling panel's lifetime cgroup peak.
- No observed contention series or newer collector-validity pass is imputed to historical runs.
- All sixteen per-cell full-pipeline costs are unavailable. Stage times are not added or substituted for them; stage peaks are neither added nor promoted to whole-run peaks.
- Original OrthoBench batch failures remain failed; separate scientific recovery does not turn those batches into clean completed runs.
- Direct source/admission/metrics pins are checked, not every transitive raw artifact or historical process-accounting assumption.
