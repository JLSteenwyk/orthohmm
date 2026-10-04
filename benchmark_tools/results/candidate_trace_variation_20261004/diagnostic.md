# Candidate Trace Variation

Retained accepted merges; no new benchmark inference or accuracy scoring.

| Run | Common merges / 8440 | Original-only merges, round 0 / 1 | Changed common anchors | All changed common anchors at cap | Max common support delta |
| ---: | ---: | --- | ---: | --- | ---: |
| 8 | 8428 | 4 / 8 | 2 | True | 2.13e-14 |
| 16 | 8439 | 0 / 1 | 1 | True | 2.13e-14 |
| 24 | 8426 | 4 / 10 | 4 | True | 2.13e-14 |

## Frozen-Engine Fixture

Nineteen genes; nine equally supported singleton satellites and a ten-gene anchor; original satellite parameters.
| Case | Unattached satellite | Merges | Rounds | Largest score perturbation |
| --- | ---: | ---: | ---: | ---: |
| baseline | [8] | 8 | 2 | 0 |
| cluster_order_only | [0] | 8 | 2 | 0 |
| one_ulp_scores_only | [7] | 8 | 2 | 2.22e-16 |

## Limits

- All accepted merge records are compared by iteration and named memberships, not only group counts or numeric labels. Rejected candidates and complete hit arrays are not replayed.
- Common accepted cluster IDs remain unchanged; the cluster-order fixture is a sensitivity result, not proof of historical relabeling.
- Capped selection with tied/near-tied evidence is consistent with the retained differences. The fixture proves this is a sufficient mechanism, not a causal explanation of every historical difference or of score-bit provenance.
- The 19-gene fixture uses the frozen source and original satellite parameters under the recorded reporting runtime, not full native benchmark execution or the original native deployment.
- Structural reference exclusion/equality is inherited from the pinned linkage report, not new benchmark scoring. Native pair files, trees and alignments remain unexamined here.
- Do not change the frozen scientific method, round scores opportunistically, choose matching repeats or claim universal invariance. A future determinism repair needs prospective numerical/ordering tests and new independent validation if assignments change.
- Existing timing observations retain unknown potentially method-dependent shared-host contention; this diagnostic neither corrects timing nor attributes variation to contention.
