# Admitted Top100 Sequence-Control Scores

Native assessment21829 completed0:0 in30:57; independent admission21830
completed0:0 in3:15. Status:`corrected_sequence_assessment_admitted`, with
accuracy_admitted=true. [Retained admission](qfo_sequence_top100_assessment_admission_21830.json)
SHA-256:`439bfb1499e6c1f5e357ec27dc0fd02279d920a6a1a77f9483d008da4849c911`.

| Endpoint | Top100 score |
| --- | ---: |
| GO similarity | 0.47979121 |
| EC similarity | 0.87497864 |
| VGNC F1 | 0.6041722529093152 |
| SwissTrees F1 | 0.6291893836042818 |
| TreeFam-A F1 | 0.5739167376813594 |
| FAS | 0.7167336149647486 |
| Project-defined secondary mean | 0.6464636398599508 |

The [updated generated table](qfo_sequence_scores_20260918_v2/scores.md)
includes both admitted sequence-control arms; its exporter rechecks the
retained admission bindings, source artifacts and native metric content.
The priorv1 table is preserved as an earlier snapshot, not overwritten.

Top100 retains11,285,357expected/emitted/mapped cross-species group-derived
clique pairs with zero mapping losses. It uses post-search truncation of
DIAMOND hits followed by frozen OrthoHMM graph/grouping code. It is not a
standalone competing tool or native phylogenetic pair predictor. Capping
hits after the search does not measure a faster search configuration.

Both arms now have validated point estimates. GO/EC/FAS are not F1, and the
secondary six-metric mean is not an official QfO aggregate. These small
point differences do not establish significance; FAS has stochastic sampling.
Neither equal E-value cutoffs nor these endpoints establish matched search
sensitivity/computational effort. The corrected HMM baseline and frozen
paired SwissTrees uncertainty analysis remain pending, so no HMM advantage
is established by this table. Shared-host job durations are workflow history,
not controlled tool-speed comparisons.
