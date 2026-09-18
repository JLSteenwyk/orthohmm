# Admitted All-Hit Sequence-Control Scores

All-hit native assessment21825 completed0:0 in31:09; independent admission
21826 completed0:0 in3:15. Status:`corrected_sequence_assessment_admitted`,
accuracy_admitted=true. [Admission report](qfo_sequence_all_hits_assessment_admission_21826.json)
SHA-256:`dc29e5410a274b5e633674bdf42c6578893bbc47912adc7b23ca6893ed917cd2`.

The[generated endpoint table](qfo_sequence_scores_20260918_v1/scores.md)
rechecks admitted hashes, conversion identity and native metric content.

| Endpoint | Score |
|---|---:|
| GO similarity | 0.47932664 |
| EC similarity | 0.87466218 |
| VGNC F1 | 0.6020310955333872 |
| SwissTrees F1 | 0.6280916644372869 |
| TreeFam-A F1 | 0.5763063958696845 |
| FAS | 0.7154481088979685 |
| Project-defined secondary six-metric mean | 0.6459776807897212 |

This control uses DIAMOND search and the frozen OrthoHMM downstream graph
and grouping workflow, with11,300,151cross-species group-derived clique
pairs and zero mapping losses. GO/EC/FAS are not F1. This is neither native
phylogenetic pair prediction nor a standalone competitor configuration.
No paired uncertainty or matched corrected HMM baseline is present yet.
Equal E-value cutoffs do not establish equivalent sensitivity or effort.

Top100conversion21828 completed0:0 in2:36, with11,285,357expected/emitted/
retained pairs and zero mapping losses. [Conversion report](qfo_sequence_top100_pairs_21828.json)
SHA-256:`955f6d151f5a276c6c4f2c4586aaffe9bfb82bb084c7bbbcb676b434fdfbd5ae`.
Its pair payload is174,201,050bytes, SHA-256
`12f78d593c7bfc0ff9b8766e01a2d391066cd579c2a5d050a04a94b3b1a2256a`.
Scoring21829 is RUNNING; admission21830 waits for it. Pending table values
are not zeros and make no scheduler-state assertion.

45focused exporter/native-validation/admission tests pass, including the
retained real score table, rejected malformed bindings, native replay
disagreement and missing/duplicate variant handling. No scientific endpoint,
parameter or inclusion rule was changed after inspecting these scores.
