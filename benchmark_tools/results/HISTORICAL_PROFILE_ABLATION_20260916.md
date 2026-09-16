# Historical Profile And Refinement Ablation

Development-exposed, descriptive controls; not the completed publication factorial.

All four partitions cover 251,378 input genes. The final replay partition
is byte-identical to the corrected fresh high-sensitivity production output.

| Stage | Groups | F1 (%) | Precision (%) | Recall (%) |
| --- | ---: | ---: | ---: | ---: |
| multipass | 51181 | 66.279051 | 73.278341 | 60.500277 |
| multipass_refined | 63245 | 69.763388 | 78.868592 | 62.542944 |
| strict_profiles | 50894 | 66.825390 | 73.657420 | 61.153180 |
| strict_profiles_refined | 62885 | 70.358998 | 78.949544 | 63.454479 |

| Descriptive effect | F1 difference (pp) | Precision difference (pp) | Recall difference (pp) |
| --- | ---: | ---: | ---: |
| profile_without_refinement | +0.546339 | +0.379079 | +0.652904 |
| profile_with_refinement | +0.595610 | +0.080952 | +0.911535 |
| refinement_without_profile | +3.484337 | +5.590251 | +2.042667 |
| refinement_with_profile | +3.533607 | +5.292124 | +2.301298 |

## Interpretation

Profile expansion has a modest positive observed effect in this historical panel.
The larger refinement effect must not be attributed to profile HMM expansion.
All arms retain the HMM-based initial search. No significance or independent
generalization claim is made from these descriptive differences.

## Limitations

- Historical development-exposed controls; descriptive effects, not prospective confirmation.
- Profile omission retains the HMM-based initial search; this is not an HMM-free comparison.
- Byte equality validates this historical endpoint, not all intermediate stages or future source versions.
- The cache hash is verified without deserializing or rechecking each normalized search hit.
- Current-source replay, full expansion/reconciliation factorial, and matched sequence-search control remain required.
- Historical source dirtiness is preserved in provenance; no clean-source equivalence is inferred.
- Replay timings are cumulative incremental costs, not independently timed ablation arms.
