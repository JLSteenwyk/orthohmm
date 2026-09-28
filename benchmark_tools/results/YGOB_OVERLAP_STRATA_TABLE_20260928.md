# YGOB Overlap-Stratified Group Recovery

Descriptive retained-count strata; not independent confirmation or subset-rescored orthology.

| Screen stratum | Pillars | Reference genes | Truth pairs | Singleton pillars |
| --- | ---: | ---: | ---: | ---: |
| screen_positive | 6952 | 77287 | 578995 | 2017 |
| screen_negative | 3298 | 6104 | 18095 | 2893 |

| Screen stratum | Method | F1 (%) | Precision (%) | Recall (%) | Gene coverage (%) | Exact pillars |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| screen_positive | OrthoHMM satellite_v2 | 93.373647 | 95.777901 | 91.087142 | 100.000000 | 3613 |
| screen_positive | OrthoHMM high sensitivity | 82.965689 | 73.209851 | 95.721379 | 100.000000 | 3199 |
| screen_positive | OrthoFinder 3.1.5 full | 95.008279 | 92.668377 | 97.469408 | 100.000000 | 3911 |
| screen_positive | OrthoFinder sequence-only checkpoint (diagnostic) | 88.625612 | 80.455002 | 98.643339 | 100.000000 | 3755 |
| screen_negative | OrthoHMM satellite_v2 | 57.932798 | 55.916909 | 60.099475 | 100.000000 | 1447 |
| screen_negative | OrthoHMM high sensitivity | 50.836356 | 46.710800 | 55.761260 | 100.000000 | 1457 |
| screen_negative | OrthoFinder 3.1.5 full | 45.016506 | 30.743257 | 84.028737 | 100.000000 | 1237 |
| screen_negative | OrthoFinder sequence-only checkpoint (diagnostic) | 37.530906 | 24.087265 | 84.935065 | 100.000000 | 1222 |

All original cross-pillar FP allocations remain, including pairs connecting the two strata.
Coverage means presence in a prediction group, not correct co-membership.
No-hit status does not establish family independence; no confidence intervals or significance claims.
All eight stratum-level ratios are defined; singleton pillar counts expose the zero-truth-pair composition.
