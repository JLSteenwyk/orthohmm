# Third Original-Release Reconciliation-On Assessment

## Admission

Cell `p1_c0_r1`, assessment job 21711, completed with exit 0:0 in 54:16.
Independent admission 21712 completed with exit 0:0 in 15 seconds.
These are assessment/validation durations, not inference timings. The
unchanged frozen admission program was run again against the same pair
manifest and terminal job; its assessment object is identical.

Retained report: `qfo_factorial_assessment_p1_c0_r1_20260918.json`.
SHA-256: `f77ca80bbabf8485fa1df10e17eda051cc8c0fbfaebd19ae12e5d698344a7940`.
This binds 4,959,440 retained pairs to all six native endpoints and execution
provenance. All 15 workflow tasks completed; FAS accounts for 52m23s of the
native trace duration. No timeout/restart or score imputation occurred.
Thirty-eight admission/native-assessment unit tests passed.

## Individual Endpoints

| Endpoint | Matched p1_c0_r0 | p1_c0_r1 |
| --- | ---: | ---: |
| GO similarity | 0.472350280 | 0.489823080 |
| EC similarity | 0.931347280 | 0.968116710 |
| VGNC harmonic P/R | 0.667694783 | 0.896753531 |
| SwissTrees harmonic P/R | 0.673850856 | 0.778297624 |
| TreeFam-A harmonic P/R | 0.576357812 | 0.566794937 |
| FAS | 0.774633359 | 0.782481533 |
| Project-defined secondary mean | 0.682705728 | 0.747044569 |

Native SwissTrees recall/precision are 0.66029942/0.94764580;
TreeFam-A recall/precision are 0.40249391/0.95775956. The TreeFam-A decrease
is preserved rather than obscured by the secondary mean. The mean is not
official QfO F1. These are point estimates, not paired significance tests;
native endpoint error fields are not substituted for paired uncertainty.

## Scope and Next Actions

This is the original-release experiment, subject to the documented Xenopus
input/reference mismatch. It must not be relabeled as corrected-input
evidence. Corrected inference remains a separate workflow.

Seven of eight original factorial cells now have admitted assessments.
The final cell `p1_c1_r1` remains in reconciliation as 21671_3, with native
admission 21673_3 and pair conversion 21675_7 waiting. After those stages
complete, submit its assessment with the actual validated pair-manifest
checksum. Run the frozen eight-cell SwissTrees counts, 42-endpoint paired
bootstrap and interaction analysis only after all eight assessments pass.
No interim tuning or multiplicity changes are authorized.

Previous turn made progress by preparing corrected comparator conversion
(`01104e0`). This turn validates and retains a newly completed empirical
assessment. Corrected comparators, full factorial uncertainty, dedicated
resource admission and publication deliverables remain open.
