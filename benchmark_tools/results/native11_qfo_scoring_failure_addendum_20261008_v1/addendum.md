# Native QfO Scoring Failure Addendum

This new evidence update supplements the frozen manuscript v4; it does not replace its bytes or four admitted-cell scores.

| Cell | Inference | Conversion | Scoring | Admitted Score | Six-Metric Mean | Submitted Pairs | Relation Coverage | Scoring CPU Slots | Scoring Memory Limit GiB |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| p1_c1_r0 | successful_composed_terminal_review | successful_unscored_group_cliques | OUT_OF_MEMORY | Unavailable | Unavailable | 11755521 | 0.5946123354776824 | 8 | 32.0 |

Inference job 23985 completed and passed its explicit full composed review. Conversion job 24035 succeeded. Assessment job 24038 ended OUT_OF_MEMORY under 8 CPUs / 32 GiB. Its FAS task exited 137; Slurm reported an OOM kill. Five other endpoint tasks completed, but the six-endpoint assessment did not complete. No score or secondary mean is admitted.

Conversion retained 11,755,521 cross-species group-clique pairs, with no mapping loss. Relation coverage is 585,180/984,137 (59.4612%). This is coverage, not accuracy, and not native phylogenetic-pair semantics.

Reviewed inference observations: 60568.563058 wall seconds, 1607563.508002 native-task CPU seconds, 20407123968 bytes native-step lifetime peak memory, including the documented wrapper/launcher scopes. These are not scoring measurements. The scoring scheduler elapsed time is 00:16:06; its exact peak RSS and the failing allocation are unavailable.

Shared Threadripper measurements while other analyses were running. CPU, memory-bandwidth and I/O contention may affect elapsed times by an unknown, potentially tool-dependent amount. These are observed shared-host timings, not estimates of isolated performance.

The failed scoring attempt is retained without retry, missing-score imputation or partial six-metric mean. A 128-GiB limit is prospectively specified for the genuinely unrun native12 scoring workflow only; it does not change inference, endpoint definitions or FAS sampling and does not guarantee success.

Evidence:
- failure: `/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmark_tools/results/native11_composed_assessment_failure_24038_20261008_v1.json` (SHA256 `b40f8075fd0168e3e7e1b895bcdba2576e2b903cc58e3484a7f3ce8b53b93e6b`)
- review: `/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/native11_composed_terminal_review_20261008_v1/review.json` (SHA256 `a87d65439a3f65dd80b199602231e8f95d10f26d2dbbee54c2f0d0eaa50abf2e`)
- conversion: `/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/results/native11_composed_qfo_pairs_v1/results.json` (SHA256 `77d33a77eb5eb1c85736ba6e7b303aa269907c8a930d74647d3b8bddff36fac8`)
- manuscript: `/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmark_tools/results/PUBLICATION_MAIN_TEXT_20261007_v4.md` (SHA256 `0b7012e3dd46bd3b028857bdb3cbe3a9e9048936b90f3df9e4833cb342d5705a`)
