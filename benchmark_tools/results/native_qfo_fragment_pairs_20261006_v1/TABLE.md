# Native Fragment-Annotation Pair Changes

Descriptive endpoint-level counts. Unflagged is not proven complete; baseline-only treats later-version annotations as missing.
No pair-IID intervals, significance, native F1 substitution or causal explanation.

| View | Endpoint bin | Relations | R0 TP | Removed TP | TP removal % | R0 FP | Removed FP | FP removal % |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| historical | annotation_positive | 527 | 42 | 0 | 0.000 | 36 | 35 | 97.222 |
| historical | all_matched_unflagged | 10238 | 3081 | 334 | 10.841 | 1672 | 1654 | 98.923 |
| historical | missing_without_positive | 0 | 0 | 0 | NA | 0 | 0 | NA |
| baseline_only | annotation_positive | 527 | 42 | 0 | 0.000 | 36 | 35 | 97.222 |
| baseline_only | all_matched_unflagged | 9681 | 2960 | 323 | 10.912 | 1538 | 1521 | 98.895 |
| baseline_only | missing_without_positive | 557 | 121 | 11 | 9.091 | 134 | 133 | 99.254 |

All 16 transition cells per view/bin, both cell marginals and the full 10,765-pair ledger are retained in report.json/pairs.tsv.
Annotation selection/history admission is inherited; selected raw entries are newly rechecked. Parsers share Bio.SwissProt.
Initial HMM search on in both native cells; failed R1 timing remains ineligible. Shared-host postprocessing costs are not inference timing; contention effects are unknown and potentially tool-dependent.
