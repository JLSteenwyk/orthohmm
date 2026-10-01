# Complete Corrected QfO Parameter Neighborhood

The frozen seven-arm analysis is complete. Job 22395 parent/batch completed
0:0 on bizon, 2 CPUs/64 GiB, elapsed 1:04. This is shared-host analysis time,
not native inference or controlled timing. The [actual submission](qfo_private_cpm_parameter_uncertainty_submission_22395.json)
binds the clean da587eed executor, four Git-verified source records, reviewed
inventory and independently admitted recovered high-CPM evidence.

The full analysis report is
`benchmarks/work/qfo_private_cpm_parameter_uncertainty_20261001.json`,
10,638,421 bytes, SHA256
`6fff1d7560c10983d8a3a2df1debf4074e1c69a88d962f4018a5de769ce9a14c`.
The [compact numerical summary](qfo_private_cpm_parameter_result_22395.json)
preserves all points, contrasts, interval/direction counts and evidence pins.
The unchanged [independent reproducer receipt](qfo_private_cpm_parameter_reproduction_20261001.json)
confirms all 18 endpoints to absolute tolerance 1e-12. It reconstructs arithmetic
from admitted counts using the same NumPy generator/quantiles, not a third
native prediction/reference audit or an independent statistical implementation.

## SwissTrees Results

F1 is recomputed from macro precision and recall in every paired bootstrap
draw, not a mean of family F1. Units below are percentages/percentage points.
The adjusted intervals retain all 18 planned endpoints, 100,000 shared family
draws, seed 20260925, linear quantiles and the same 18 reference families.

| Arm | F1 (%) | F1 difference (pp) | Adjusted interval (pp) | Family wins/ties/losses |
| --- | ---: | ---: | --- | --- |
| Frozen control | 83.351322 | Reference | Reference | Reference |
| CPM 0.08 | 79.101373 | -4.249949 | [-21.939484, 0.788718] | 3/13/2 |
| CPM 0.12 | 83.253016 | -0.098306 | [-1.050447, 0.747525] | 1/15/2 |
| Minimum normalization 0.024 | 83.351322 | 0 | [0, 0] | 0/18/0 |
| Minimum normalization 0.036 | 83.351322 | 0 | [0, 0] | 0/18/0 |
| Minimum margin 1.2 | 82.839072 | -0.512250 | [-4.117304, 1.352279] | 1/16/1 |
| Minimum margin 1.8 | 83.564753 | 0.213431 | [-0.928477, 1.957055] | 1/16/1 |

High-CPM precision/recall are 95.460954%/73.813450%. Their differences are
-0.056721/-0.120676 pp, adjusted intervals [-0.264550, 0.007904] and
[-1.501235, 1.107420]. ASTER and CASP lose family F1; BAR improves and the
other 15 tie. This is descriptive localization, not a causal explanation.
All 18 adjusted intervals include zero. No statistically supported parameter
improvement is established, and intervals including zero do not prove equivalence.
Zero observed normalization contrasts do not establish identical whole-proteome
predictions or equivalence on unobserved families. No defaults are changed.

The other six point estimates and five contrasts exactly match the retained
partial report. The formerly missing high arm now has separate private-runtime
recovery, native-pair conversion and [independent score admission](QFO_PRIVATE_CPM_SCORE_RESULT_22394.md).
Historical SIGSEGV/allocator/admission failures remain failed; recovery does not
explain their cause or prove memory safety. A content-equivalent private control
was validated before recovered inference. This is not a controlled runtime study.

## Artifacts And Limits

[Complete native-score table](qfo_parameter_complete_export_20261001/scores.md),
[all 18 intervals](qfo_parameter_complete_export_20261001/intervals.tsv), and
[SwissTrees figure](qfo_parameter_complete_export_20261001/qfo_parameter_neighborhood.pdf)
are generated from the exact admitted result/reproduction pair. Full manifests
remain locally retained; compact readback receipts bind their bytes and outputs.
The earlier [partial figure](qfo_parameter_partial_export_20261001_v2/qfo_parameter_neighborhood.pdf)
preserves its then-missing high arm and explicit historical-source amendment.

QfO is development-exposed; the panel is not independent confirmation, an
OrthoFinder superiority test or broad evolutionary robustness. Only SwissTrees
has family-level intervals here. GO/EC/FAS are not F1, the six-score mean is
secondary, and unseeded FAS/dependence/attrition remain. Other QfO uncertainty,
controlled scaling, rights and the final publication package remain incomplete.
