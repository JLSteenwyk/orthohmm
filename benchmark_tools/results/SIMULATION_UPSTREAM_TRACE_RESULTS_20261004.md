# Retained Simulation Upstream Trace

Exploratory finite-panel counts, not causal effects or independent confirmation.
All seven conditions include all ten fixed seeds. Direct significant hits, graph edges
and graph paths are not native ortholog predictions. Native TP/FN and cross-candidate
counts reproduce the preceding frozen scores/oracle readback.

| Condition | True pairs | Across candidates | Different graph components | Connected but separated | Direct edge but separated | Native FN within candidates |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| baseline | 29,051 | 177 | 177 | 0 | 0 | 98 |
| divergent | 29,051 | 12,736 | 11,442 | 1,294 | 117 | 60 |
| divergent_turnover | 30,948 | 13,072 | 11,554 | 1,518 | 128 | 142 |
| missing20 | 18,652 | 138 | 138 | 0 | 0 | 60 |
| taxon_count_control | 12,458 | 92 | 92 | 0 | 0 | 43 |
| turnover | 30,948 | 196 | 196 | 0 | 0 | 206 |
| uneven_taxa | 12,419 | 162 | 162 | 0 | 0 | 43 |

## Cross-Candidate Hit Evidence

Each cell below sums true pairs over ten seeds; hit orientations refer to significant
initial sequence-HMM hits, not every scored prefilter candidate.

| Condition | Different components: 0 / 1 / 2 orientations | Connected but separated: 0 / 1 / 2 orientations |
| --- | ---: | ---: |
| baseline | 177 / 0 / 0 | 0 / 0 / 0 |
| divergent | 11,441 / 0 / 1 | 1,177 / 1 / 116 |
| divergent_turnover | 11,553 / 1 / 0 | 1,384 / 0 / 134 |
| missing20 | 138 / 0 / 0 | 0 / 0 / 0 |
| taxon_count_control | 92 / 0 / 0 | 0 / 0 / 0 |
| turnover | 196 / 0 / 0 | 0 / 0 / 0 |
| uneven_taxa | 162 / 0 / 0 | 0 / 0 / 0 |

## Evidence Limits

- Development-exposed descriptive trace; no population CI, significance or causal intervention.
- Absent initial significant hits do not distinguish prefilter rejection, score rejection, caps or ranking.
- Saved final graph precedes refinement/candidate expansion; paths may traverse non-reference-family genes.
- Connected but separated localizes to grouping/refinement/expansion as a combined boundary, not one isolated algorithm.
- Gene-to-seed membership inside merged candidates is not retained; no pre-expansion partition is invented.
- Native counts reproduce, but graph closure and direct hits are not orthology predictions.
- No search/inference/timing rerun, new default or full publication completion.
