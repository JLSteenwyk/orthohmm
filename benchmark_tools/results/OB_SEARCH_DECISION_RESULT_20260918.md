# Observed OrthoBench Search Rejections

Frozen diagnostic 21856 completed 0:0 in 6:59; independent raw-output recount
21857 completed 0:0 in 6 seconds. All 144 species directions and 81,466
distinct directed within-reference-family pairs were checked.

| Observed decision | Directed pairs | Historical hit presence |
| --- | ---: | --- |
| Accepted | 31,479 | Present for all |
| Not selected by prefilter | 48,368 | Absent for all |
| Scored but E-value not below 1e-4 | 1,619 | Absent for all |
| Total | 81,466 | Zero presence disagreements |

Thus, the observed CPU diagnostic reproduces historical accepted-hit presence
on this watched set. Most observed absent hits were excluded before scoring.
This does not authenticate the original executable/runtime, establish
numerical score equality, or extend to unwatched pairs. It is not evidence
that every excluded pair should be accepted: the descriptive pair universe
includes within-species and low-certainty members and is not the official
weighted recall statistic. No counterfactual scoring of excluded pairs was
performed. Graph paths and later processing can group genes without a direct
hit, so rejection-stage localization is not a causal explanation of final F1.

## Evidence

- [Frozen protocol](OB_SEARCH_DECISION_PROTOCOL_20260918.md).
- [Independent recount](ob_search_decisions_audit_20260918.json), unmodified
  from `benchmarks/work/ob_search_decisions_audit_20260918.json`, 88,615 bytes;
  SHA-256 `8081c8ff005c911b0122ba3a797aaee51c1ac65a1c3f17f551512d3cb47ce730`.
- Driver report `benchmarks/results/ob_search_decisions_v1/report.json`;
  SHA-256 `09ab10be8fba1b31917b2903836c55f70de5680e390ea2f47ab88f80a0fb8398`.
- All 144 raw NPZ candidate archives and watched-pair TSVs remain beside
  the driver report, with their hashes bound by both reports.

The independent audit reconstructs decisions from raw candidate arrays,
checks input ID order/universes, preserves reference-family overlap labels,
compares each TSV row and counter, and rehashes the evidence. It does not
independently recompute HMM scores. Neither report admits accuracy or
publication readiness. Shared-host diagnostic duration is not a comparative
end-to-end resource result. No frozen benchmark prediction was replaced.
