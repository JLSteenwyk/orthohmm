# Installed OrthoBench: Search Agreement, Candidate Divergence

## Retained Evidence

The [fresh run readback](INSTALLED_ORTHOBENCH_READBACK_20260926.md) found
73.821569% versus historical 74.106074% F1. This post-hoc diagnostic compares
saved evidence without another inference run or parameter selection.

The historical, previously admitted SHA-256-pinned local pickle contains
normalized directed nonself hits. The fresh installed pipeline persists
normalized indexed scores before graph construction. Map integer indices to
gene IDs, exclude self-hits and collapse duplicate directed pairs by their
maximum score. Compare gene universes and all hits, not only reference-family
members. Do not divide the fresh scores by protein length again.

| Search comparison | Count |
|---|---:|
| Fresh checkpoint rows | 18,486,508 |
| Fresh self-hit rows excluded | 251,135 |
| Unique directed nonself hits, each run | 18,235,373 |
| Historical-only / fresh-only hits | 0 / 0 |
| Shared scores exactly equal | 18,225,151 |
| Shared scores unequal | 10,222 |
| Outside absolute/relative tolerance of 1e-12 | 0 |

Maximum absolute difference is 8.881784197001252e-16; maximum relative
difference is 4.0987776641924614e-16. These are observed numerical differences,
not evidence identifying a particular dependency or floating-point operation.
The [machine-readable result](installed_ob_search_comparison_20260926.json)
preserves a concrete differing-score witness and all input hashes.

## Where Agreement Ends

The candidate-family partitions, before phylogenetic inference, differ:

| Candidate partition comparison | Count |
|---|---:|
| Historical groups | 54,445 |
| Fresh groups | 54,495 |
| Identical groups, ignoring labels | 51,872 |
| Historical-only groups | 2,573 |
| Fresh-only groups | 2,623 |
| Genes in changed groups | 30,881 |

Thus the full-run disagreement cannot be attributed solely to phylogenetic
reconciliation, and there is no observed accepted-hit presence loss between
these retained search artifacts. The next unresolved interval is graph
construction, clustering, profile expansion, refinement and candidate merging.
Tiny score changes can potentially affect ties or thresholds, but their causal
contribution has **not** been tested. Different graph/clustering dependencies,
ordering and execution paths remain possibilities, not established causes.

The fresh run overwrites early cluster files as it proceeds, so saved final
candidate outputs alone cannot identify the first differing clustering step.
Next compare derived graph edges under the same implementation using both
retained hit sets, before considering any explicitly scoped native diagnostic.
Do not launch another complete run merely to seek a matching score.

## Reproduction and Limits

```sh
/home/bizon/anaconda3/bin/python -m benchmark_tools.compare_installed_ob_search \
  --repo . --output /tmp/installed_ob_search_comparison.json
/home/bizon/anaconda3/bin/python -m pytest -q \
  tests/unit/test_compare_installed_ob_search.py \
  tests/unit/test_audit_installed_orthobench.py
```

All 43 focused tests pass, including 17 new tests for directed identity,
duplicate maxima, self-hit exclusion, numerical differences, malformed arrays
and strict partition coverage. The reader refuses an existing output and
rechecks input hashes after comparison. Deserialization is restricted to the
already admitted historical cache hash; do not substitute an untrusted pickle.

Raw candidates, rejected scores and E-values are not compared. Equal retained
hit presence is not proof of raw search equivalence or biological accuracy.
No historical score, default, scientific executor, or benchmark table changed.
