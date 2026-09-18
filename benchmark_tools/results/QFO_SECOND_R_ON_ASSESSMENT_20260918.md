# Second Original-Release Reconciliation-On Assessment

## Admission

Cell `p0_c1_r1`, assessment job 21703, completed successfully in 36:30;
validator 21704 completed successfully in 15 seconds. These elapsed
times cover scoring and validation, not method inference.
The frozen independent admission program was rerun and produced an
identical assessment object. Admission artifact:
`qfo_factorial_assessment_p0_c1_r1_20260918.json`, SHA-256
`dd5236702a576dc278787c0d94385b0a426428ba383392a8d051e43bab46ce27`.
It binds the 5,580,807 retained pairs to native scoring outputs and
execution provenance. It is not a corrected-release result.

## Individual Endpoints

| Endpoint | Matched p0_c1_r0 | p0_c1_r1 |
| --- | ---: | ---: |
| GO similarity | 0.471009680 | 0.490377300 |
| EC similarity | 0.919343710 | 0.969003540 |
| VGNC harmonic P/R | 0.647990059 | 0.900300866 |
| SwissTrees harmonic P/R | 0.676019536 | 0.807216849 |
| TreeFam-A harmonic P/R | 0.604915207 | 0.580035916 |
| FAS | 0.746946524 | 0.761326827 |
| Project-defined secondary mean | 0.677704120 | 0.751376883 |

SwissTrees native recall/precision are 0.70006819/0.95309196;
TreeFam-A recall/precision are 0.41539981/0.96085201. The TreeFam-A
decrease is retained alongside the other point-estimate increases.
The secondary mean is not official QfO F1. No paired confidence interval
or superiority conclusion is implied by this interim table. Native
endpoint error fields have endpoint-specific semantics and must not be
substituted for paired method-difference uncertainty.

## Next Actions

Finish the remaining two reconciliation-on cells and their downstream
admissions, then run the prespecified complete eight-cell SwissTrees
family bootstrap and interaction contrasts. Do not select parameters or
alter multiplicity handling from these interim observations. All values
here remain subject to the documented original-release Xenopus mismatch;
fresh corrected-input inference has separately started as array 21706.
