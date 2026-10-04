# Selected QfO OrthoHMM Stage Provenance

All 24 earlier rows and their scores remain unchanged. New records are incremental, not full inference.
The same checked replay observation appears in both rows; it is not two independent runs.

| Method | Cell | Cached Replay Worker (s) | Arm Preparation (s) | Reconciliation (s) | Full Pipeline |
| --- | --- | ---: | ---: | ---: | --- |
| OrthoHMM high sensitivity | p1_c0_r0 | 2987.810 | 3.036 | NA (R off) | NA |
| OrthoHMM phylogeny satellite_v2 | p1_c1_r1 | 2987.810 | 77.641 | 6822.365 | NA |

## Limits

- Two selected QfO OrthoHMM chains are linked through exact score/conversion/candidate/replay/admission metadata. Raw predictions, FASTAs, numeric hits and source/runtime inventories are inherited from prior admissions, not re-audited here.
- Five stage associations add a shared replay worker in both rows, two shared candidate-arm intervals and one reconciliation interval. These are not five new runs or independent repetitions. The replay wall/CPU is a checked worker wrapper; internal replay timings are a distinct scope.
- Initial HMM search is absent from cached replay costs. No stage sum, memory sum or full-pipeline cost is inferred. Full costs for the exact selected cached executions remain unavailable; upstream native equivalence is not used as a cost substitute.
- Candidate preparation CPU and memory, and native conversion costs are unestablished. Recorded argv is not a new executable/dependency or historical consumption attestation. All shared-host distortion remains unknown and potentially method dependent.
