# YGOB Curated Group Recovery

Completion gates require separate verification; this report alone does not establish publication readiness.

| Method | F1 (%) | Precision (%) | Recall (%) | Reference-gene coverage (%) | Predicted-group coverage (%) | Exact reference groups |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| OrthoHMM satellite_v2 | 92.233654 | 94.418045 | 90.148051 | 100.000000 | 99.990078 | 5060/10250 |
| OrthoHMM high sensitivity | 82.038608 | 72.474709 | 94.510375 | 100.000000 | 100.000000 | 4656/10250 |
| OrthoFinder 3.1.5 full | 92.318524 | 88.017010 | 97.062084 | 100.000000 | 100.000000 | 5148/10250 |
| OrthoFinder sequence-only checkpoint (diagnostic) | 85.572844 | 75.806425 | 98.227905 | 100.000000 | 100.000000 | 4977/10250 |

Predicted-group coverage is the fraction of supplied groups retaining at least one scored gene.

## Paired Differences Versus Full OrthoFinder

20,000 paired pillar replicates; PCG64 seed 20260917. Differences and intervals are percentage points.

| Contrast | Metric | Difference | Nominal 95% CI | Bonferroni CI (six contrasts/metrics) |
| --- | --- | ---: | --- | --- |
| Primary: OrthoHMM satellite_v2 | f1 | -0.084870 | [-0.487503, 0.312676] | [-0.622528, 0.445222] |
| Primary: OrthoHMM satellite_v2 | precision | 6.401035 | [5.983125, 6.839557] | [5.841375, 6.993930] |
| Primary: OrthoHMM satellite_v2 | recall | -6.914033 | [-7.542547, -6.300239] | [-7.759841, -6.091851] |
| Secondary: OrthoHMM high sensitivity | f1 | -10.279916 | [-11.021330, -9.540287] | [-11.260924, -9.298064] |
| Secondary: OrthoHMM high sensitivity | precision | -15.542301 | [-16.695003, -14.370327] | [-17.079385, -13.991748] |
| Secondary: OrthoHMM high sensitivity | recall | -2.551709 | [-3.000514, -2.098126] | [-3.151369, -1.956990] |

## Limitations

- Curated homolog-group co-membership, including within-species pairs, not resolved pairwise orthology.
- Novel-taxon transfer does not establish family-disjoint or unrestricted generalization.
- Predictions outside the reference universe are projected out, not counted as false positives.
- Coverage and missing predictions must be interpreted alongside F1.
- Run completion, frozen inputs and versions, overlap and resource audits require separate verification.
- Approximate intervals assume exchangeable reference pillars.
- Cross-pillar false positives are allocated half to each endpoint pillar.
- Group recovery does not establish resolved pairwise orthology.
