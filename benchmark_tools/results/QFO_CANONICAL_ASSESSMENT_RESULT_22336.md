# Canonical QfO Assessment Result

## Completed And Independently Validated

Scoring job 22336 completed 0:0 in 00:30:00 with eight allocated CPUs.
The independent admission validator checked the exact plan/submission,
conversion, frozen environment, source/input identities, execution receipts,
complete native output inventory and all 15 fresh successful workflow tasks.
All six native endpoint outputs and SwissTrees family assessments passed
content and aggregation validation. Full local admission:
`benchmarks/work/qfo_canonical_assessment_20260927/admission.json`,
577,996 bytes, SHA256
`194d38a5ef6448e9f2b81a1def944a51dd16885037f73e180b0330b4feb51cd6`.

The [historical comparator was independently readmitted](qfo_canonical_historical_readmission_20260927.json).
The [generated comparison](qfo_canonical_comparison_20260927/scores.md) and
[machine-readable results](qfo_canonical_comparison_20260927/results.json)
retain both arms, native axes, available precision/recall, assessed-relation
counts, admission identities and the exporting source identity.

## Interpretation

GO, EC, VGNC, SwissTrees and TreeFam-A endpoint scores are exactly unchanged.
FAS is 0.762588129455 versus historical 0.762993312100, a difference of
-0.000405182646. The native FAS scorer uses unseeded sampling: this difference
cannot be attributed solely to ordering and is not evidence of an ordering
accuracy regression. No preferred sample was selected or scorer retried.

The secondary six-metric mean is 0.761442488335 versus 0.761510018776.
This is a project-defined mean of heterogeneous endpoints, not F1 or an
official universal ranking; its change is entirely from FAS. Native uncertainty
fields are not new paired confidence intervals. The 5,959,535 submitted
canonical pairs had zero mapping loss; native NR_ORTHOLOGS axes are
challenge-specific and must not all be interpreted as total pair coverage.

These results complete the frozen downstream ordering assessment, not fresh
end-to-end QfO search, independent generalization, controlled timing or all
publication requirements. Historical scores remain retained as their own
evaluated predictions. Production defaults are unchanged. The prior
[prediction differences](QFO_CANONICAL_RESULT_22333.md) remain real despite
unchanged scores on five endpoints. Do not choose ordering based on scores.
