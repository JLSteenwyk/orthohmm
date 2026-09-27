# Matched-Graph Stage Trace Results

All 70 cells were traced under the
[post-hoc protocol](MATCHED_GRAPH_STAGE_TRACE_PROTOCOL_20260926.md), committed
as `3915615` before stage outcomes were computed. The implementation was
committed as `2f964fd` before execution. This localizes the observed difference;
it does not isolate a causal mechanism or introduce confirmatory endpoints.

## Findings

Equal-dataset mean pair F1 (%) at the partition stages:

| Stage | HMM | DIAMOND |
|---|---:|---:|
| Initial clustering | 83.9344 | 80.6939 |
| Singleton-augmented clustering | 83.9571 | 80.7386 |
| Final refinement | 83.6097 | 80.5180 |

Most of the final between-arm difference is already present at initial
clustering. In the two divergent conditions, initial-clustering F1 is
59.0573/49.3443 and 59.9822/48.8569 (HMM/DIAMOND). Subsequent singleton
attachment does not account for most of these differences.

Refinement reduces the overall mean F1 in both arms. Across the 35 dataset
instances it adds 26 true and 1,221 false pairs for HMM, and 8 true and 895
false pairs for DIAMOND, removing none. These are descriptive pair counts,
not independent observations or a reason to retune the frozen method on this
panel. The negative result is retained.

The final HMM partitions contain a mean 170.54 correct pairs per dataset
without a direct accepted hit in either direction, versus 50.20 for DIAMOND.
They also contain more false pairs without a direct hit (48.74 versus 26.31).
Indirect connections can support assignments; absence of direct evidence is
not proof of a search failure or unsupported graph assignment.

The full table also reports direct-hit, graph-edge and component-closure
overlap with orthology truth. Those objects are homology/graph evidence, not
native ortholog predictions, so their overlap F1 must not be relabeled as
tool accuracy. Different hit identities, weights and rankings remain coupled.

## Verification and Reproduction

- [All stages and conditions](matched_graph_stage_trace_v2_20260926/results.md)
- [Every cell, transition count, support count and evidence checksum](matched_graph_stage_trace_v2_20260926/results.json)
- Every final-stage score exactly reproduces the frozen score artifact.
- The focused trace, audit, scoring, plotting and statistical reproduction
  suite passed 65 tests, including within-species bridges, isolates, duplicate
  directions, transition accounting, incomplete panels and altered final scores.

The first local trace remains retained but unselected. The selected v2 repeats
the same calculation after the scorer accepted the exact strengthened-readback
hash; its cell records and summaries are identical. This refreshes a source
dependency fingerprint, not scientific evidence. A full native-artifact rescore
using v2 also reproduced all original cell records, 24 contrasts and bootstrap
configuration exactly. The reproduction guide now selects v2; the old pinned
auditor source predates the current checkout and correctly fails its source
identity gate there. No historical report was overwritten.

```bash
python -m benchmark_tools.trace_matched_graph_stages \
  --readback benchmark_tools/results/matched_graph_readback_v2_20260926.json \
  --scores benchmark_tools/results/matched_graph_scores_20260926/results.json \
  --protocol benchmark_tools/results/MATCHED_GRAPH_STAGE_TRACE_PROTOCOL_20260926.md \
  --output /absolute/fresh/matched_graph_stage_trace
```

This requires retained local native artifacts. It adds neither independent
validation nor controlled resource evidence. No native inference was rerun,
defaults changed, or new confidence intervals calculated. Full-pipeline
OrthoFinder comparisons and the original claim boundaries remain unchanged.
