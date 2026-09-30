# GO/EC Pair-Composition Figure

The [generated figure](figures_qfo_scored_pair_decomposition_20260930/qfo_pair_decomposition.pdf)
and [complete table](figures_qfo_scored_pair_decomposition_20260930/scores.md)
display all seven method contrasts against full OrthoFinder 3.1.5 for GO and
EC. They reuse the [historically bound 56-comparison panel](QFO_SCORED_PAIR_TRANSITIVE_BINDING_20260927.md),
not a new inference, annotation-score calculation or uncertainty analysis.

The exporter checks all 56 comparison inventories and their exact integer
score-sum arithmetic before selecting the 14 baseline contrasts. Reversing
a stored comparison also reverses its denominators and exclusive terms.
It verifies all 16 endpoint counts/means and rechecks 54 source/input record
entries before and after export. These are direct recorded bindings, not a
fresh audit of every transitive raw annotation or prediction input.

For method A versus reference B, the rounded-mean difference is exactly
shared score sum / A's scored count minus shared score sum / B's scored count,
plus A-only score sum / A's count, minus B-only score sum / B's count.
The validation uses rational arithmetic; display values are floating point.
Shared-pair scores are identical at the retained six-decimal precision in
all 56 comparisons. That does not prove equality of unavailable full-precision
values or annotation/scorer correctness.

| Phylogenetic OrthoHMM minus full OrthoFinder | GO | EC |
| --- | ---: | ---: |
| OrthoHMM scored pairs | 84,211 | 117,460 |
| OrthoFinder scored pairs | 163,557 | 175,361 |
| Shared pairs | 75,894 | 101,712 |
| Shared fraction of OrthoHMM pairs | 90.12% | 86.59% |
| Shared fraction of OrthoFinder pairs | 46.40% | 58.00% |
| Shared denominator term, score points | +21.571332 | +27.772815 |
| OrthoHMM-only term, score points | +4.569645 | +12.451329 |
| Negative OrthoFinder-only term, score points | -24.060846 | -37.272115 |
| Net difference, score points | +2.080130 | +2.952029 |

One score point is 0.01 similarity units; these are not F1 differences. The
relatively small net differences coexist with much larger opposing terms.
This is an exact eligible-pair composition identity at serialized precision,
not causal attribution, a selection-bias diagnosis or a superiority test.
Scored-pair fractions are not proteome-wide prediction coverage. Conditioning
on shared pairs changes the endpoint and is not a replacement benchmark.
The sequence-only comparator remains a pre-phylogeny MCL checkpoint;
FastOMA uses a supplied OrthoFinder species tree.

## Reproduce

```bash
python -m benchmark_tools.plot_qfo_scored_pair_decomposition \
  --results benchmark_tools/results/qfo_scored_pair_panel_20260927_v2.json \
  --sha256 7a0e281fc72143f50fbc0c419127fd8bd16c12ad8cddcc57b86f09b309e377a6 \
  --output /fresh/path/qfo-scored-pair-figure
python -m pytest -q tests/unit/test_plot_qfo_scored_pair_decomposition.py \
  tests/unit/test_compare_qfo_scored_pairs.py \
  tests/unit/test_run_qfo_scored_pair_panel.py
```

43 focused tests pass. The [export manifest](figures_qfo_scored_pair_decomposition_20260930/manifest.json)
pins the five exports, input panel and checked records. The separate
[readback/visual-review receipt](qfo_pair_decomposition_review_20260930.json)
verifies all exported component sums and a readable, nonblank PNG with both
panels, labels and boundary notes. The PDF page has not had a separate visual
review. No native score, method default, timing run or historical receipt was
changed. Appropriate QfO comparison uncertainty remains unresolved.
