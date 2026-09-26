# Matched-Recall Graph Result

All 70 graph arms completed and passed independent numeric, graph-edge and
partition readback. Scoring followed the
[frozen protocol](MATCHED_SEARCH_GRAPH_PROTOCOL_20260926.md), with no changes
to thresholds, graph settings or endpoints after inspecting outcomes.

The HMM arm has mean cross-species pair F1 **83.6097%**, versus **80.5180%**
for the DIAMOND arm: **+3.0917 percentage points**. The prespecified paired
seed-block interval, adjusted for eight F1 contrasts, is **[+1.6782, +4.5313]**
points. All five overall seed effects are positive. The strongest differences
occur under divergence (+9.6223 points) and divergence with turnover (+10.6832
points). Uneven sampling slightly favors DIAMOND (-0.1014 points), with an
adjusted interval spanning zero. The other four condition intervals touch zero.
These approximate intervals have only five independent seed blocks.

The overall trade-off is mean precision **90.7136% vs 91.0967%** and recall
**83.4079% vs 79.9802%** (HMM vs DIAMOND). The exploratory marginal interval
for the precision difference spans zero. All genes are represented in both
partitions; mean nonsingleton coverage is **94.7839% vs 93.4863%**. Scores are
means of dataset statistics, not pooled-pair statistics or harmonic means of
the displayed mean precision and recall.

## Artifacts and Verification

- [Generated complete table](matched_graph_scores_20260926/results.md)
- [Machine-readable scores, per-seed effects and intervals](matched_graph_scores_20260926/results.json)
- [Figure PNG](matched_graph_figure_v2_20260926/matched_graph.png),
  [PDF](matched_graph_figure_v2_20260926/matched_graph.pdf),
  [SVG](matched_graph_figure_v2_20260926/matched_graph.svg)
- [Figure provenance](matched_graph_figure_v2_20260926/manifest.json)
- [Native readback](matched_graph_readback_20260926.json)

Scoring code was committed as `ed7f7ed` before execution. The result JSON has
SHA-256 `d31dbc0f4e3f02dfb534c8b59b559db00260227f23f2640b4fed13d5aee8ec89`.
The scoring/search/graph suite passed 110 tests; two additional figure tests
verify aligned rows and rejection of altered effects. Figure values are
recomputed from paired score records. The final v2 PNG was visually reviewed:
all labels, intervals, legends and limitations are visible without overlap.
The initial local rendering is retained but not selected for the publication
figure because its two panels had slightly different automatic row limits.
The generated manifest precedes this visual check and keeps its initial
`visual_review_complete=false`; this note records the subsequent review.

## Claim Boundary

The defensible finding is an advantage for HMM-derived initial-search evidence
in this fixed graph procedure on these matched-recall simulations, especially
the two divergent conditions. The design does not isolate a causal mechanism:
score distributions, rankings and hit identities still differ between arms.
Do not describe the gain as proof that HMMs alone cause better orthology inference.

This is not independent validation, full-pipeline superiority, or a comparison
with OrthoFinder. The simulations were development-exposed; profile expansion,
candidate expansion and phylogeny were off. The selected DIAMOND threshold
must not be transferred to real benchmarks as a claim of matched sensitivity.
Shared-host incremental resource records do not establish matched efficiency.
Scientific defaults and prior real-data results remain unchanged.

Next integrate the bounded result and figure into the manuscript and
claim-to-evidence checklist, while retaining dedicated timing, broader
uncertainty/generalization, installed full-run readback and release obligations.
