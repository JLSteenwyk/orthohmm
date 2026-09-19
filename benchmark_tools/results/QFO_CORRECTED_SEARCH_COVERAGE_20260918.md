# Corrected QfO Search Coverage

Job21793 completed0:0 in14:59 with2allocated CPUs and192GiB requested RAM
on the shared host. Frozen executor
`3a83f9db2981c7bad3d984738c599cf0d1bb4585` remains clean. The completed
source report is retained at
`benchmarks/work/qfo_corrected_hit_coverage_20260918/report.json`, SHA-256
`0e2b668fc5b2829f0a73175245e6c491986b0f2b1b89b2a4111776917d521e0f`.

This is a label-free search diagnostic, **not orthology accuracy**. The HMM
checkpoint is written before multi-sequence profile expansion in the frozen
scientific core. Its hits do not include all later profile-based recoveries,
candidate-family expansion or reconciliation predictions.

## Results

All three searches cover the same984137-protein,78-species input universe.
The [machine-readable table](qfo_hit_coverage_20260918/coverage.tsv) and
[source-bound summary](qfo_hit_coverage_20260918/summary.json) retain counts,
overlap statistics and provenance.

| Search checkpoint | Directed hits | Non-self hits | Cross-species hits | Queries without non-self hits | Queries without cross-species hits |
| --- | ---: | ---: | ---: | ---: | ---: |
| HMM initial search | 90,687,327 | 89,703,492 | 80,995,173 | 154,322 | 251,557 |
| DIAMOND all-hits | 593,510,904 | 592,530,075 | 536,375,270 | 117,956 | 228,520 |
| DIAMOND top100 | 321,164,891 | 320,184,088 | 302,314,946 | 117,956 | 228,520 |

HMM and DIAMOND all-hits share48,603,572non-self directed hits:
54.1825% of the HMM set and8.2027% of the DIAMOND set. There are41,099,920
HMM-only and543,926,503DIAMOND-only non-self hits; Jaccard overlap is0.076707.
Thus the HMM hit set is not simply a smaller subset of DIAMOND's set.

The top100 control retains54.0368% of DIAMOND's non-self hits, with zero
hits absent from the all-hits source. HMM/top100 non-self overlap is40,126,348
hits, or44.7322% of HMM hits. Top100 is a post-search per-query/per-target-
species reporting cap, not a global100-hit limit or a reimplementation of
the HMM prefilter. It preserves query-level non-self/cross-species coverage
in this dataset while dropping many hits; that does not establish preserved
biological sensitivity.

## Validation And Limits

The completed source job checks bound input evidence, canonical checkpoint
integrity and exact directed-hit intersections. The compact exporter pins
that full report's checksum, rechecks source/helper/manifest hashes, and
checks count conservation, species-direction totals, score histogram totals,
overlap arithmetic and the top100 subset relation. Its48focused tests pass:

```sh
python -m pytest -q tests/unit/test_export_qfo_hit_coverage.py tests/unit/test_summarize_sorted_search_hits.py
python benchmark_tools/export_qfo_hit_coverage.py \
  --source benchmarks/work/qfo_corrected_hit_coverage_20260918/report.json \
  --output benchmark_tools/results/qfo_hit_coverage_20260918
```

The export command requires a new output directory. Its checks do not
independently reproduce intersections from every raw hit; the full report
and arrays remain retained for that purpose. No benchmark reference labels
are used. Equal E-value cutoffs do not equate calibration, sensitivity or
computational effort. Missing hits do not distinguish prefilter rejection
from score rejection. More hits and broader query coverage are not evidence
of better precision, recall or F1. Corrected HMM downstream scoring and
paired comparisons remain pending; no method advantage follows from this
diagnostic. Shared-host diagnostic elapsed time is not inference timing.
