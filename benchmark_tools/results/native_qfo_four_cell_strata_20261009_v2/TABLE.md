# Four-Cell Native SwissTrees Fixed Strata

Scores in percent; differences in percentage points. NA is an empty bin, not zero.
R at P0/C0, C at P0/R0, P at C0/R1. Initial HMM search stays on; no subgroup intervals.

## Sequence

| Bin | Families | P0/C0/R0 F1 | P0/C0/R1 F1 | P0/C1/R0 F1 | P1/C0/R1 F1 | R F1 | R PPV | R TPR | C F1 | C PPV | C TPR | P F1 | P PPV | P TPR |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| all | 18 | 68.918 | 78.957 | 68.571 | 78.637 | +10.039 | +30.521 | -6.533 | -0.347 | -3.523 | +4.375 | -0.320 | -0.010 | -0.463 |
| concentrated | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| explicit_fragment | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| higher_entropy | 9 | 67.386 | 87.360 | 66.650 | 86.597 | +19.974 | +37.818 | -1.749 | -0.736 | -2.345 | +3.079 | -0.763 | -0.529 | -0.926 |
| lower_entropy | 9 | 68.521 | 68.941 | 69.120 | 69.074 | +0.420 | +23.224 | -11.318 | +0.599 | -4.702 | +5.671 | +0.133 | +0.509 | +0.000 |
| missing | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| missing_entropy | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| no_explicit_fragment | 18 | 68.918 | 78.957 | 68.571 | 78.637 | +10.039 | +30.521 | -6.533 | -0.347 | -3.523 | +4.375 | -0.320 | -0.010 | -0.463 |
| not_concentrated | 18 | 68.918 | 78.957 | 68.571 | 78.637 | +10.039 | +30.521 | -6.533 | -0.347 | -3.523 | +4.375 | -0.320 | -0.010 | -0.463 |
| not_short_relative | 11 | 69.008 | 80.389 | 68.772 | 79.890 | +11.381 | +29.920 | -5.121 | -0.236 | -3.961 | +5.889 | -0.499 | -0.016 | -0.758 |
| short_relative | 7 | 68.679 | 76.410 | 67.976 | 76.410 | +7.731 | +31.465 | -8.753 | -0.703 | -2.836 | +1.996 | +0.000 | +0.000 | +0.000 |

## Domain

| Bin | Families | P0/C0/R0 F1 | P0/C0/R1 F1 | P0/C1/R0 F1 | P1/C0/R1 F1 | R F1 | R PPV | R TPR | C F1 | C PPV | C TPR | P F1 | P PPV | P TPR |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| all | 18 | 68.918 | 78.957 | 68.571 | 78.637 | +10.039 | +30.521 | -6.533 | -0.347 | -3.523 | +4.375 | -0.320 | -0.010 | -0.463 |
| median_pfam_types_at_least_two | 6 | 63.986 | 84.932 | 64.629 | 84.932 | +20.946 | +42.987 | -1.262 | +0.643 | -0.958 | +3.839 | +0.000 | +0.000 | +0.000 |
| median_pfam_types_below_two | 12 | 71.006 | 75.897 | 70.334 | 75.401 | +4.891 | +24.288 | -9.169 | -0.672 | -4.806 | +4.643 | -0.496 | -0.015 | -0.694 |
| repeated_type_fraction_at_least_quarter | 3 | 65.595 | 74.731 | 66.582 | 74.731 | +9.136 | +37.390 | -13.296 | +0.987 | -1.054 | +4.496 | +0.000 | +0.000 | +0.000 |
| repeated_type_fraction_below_quarter | 15 | 69.545 | 79.749 | 68.954 | 79.371 | +10.204 | +29.147 | -5.181 | -0.591 | -4.017 | +4.351 | -0.378 | -0.012 | -0.556 |

## Duplication

| Bin | Families | P0/C0/R0 F1 | P0/C0/R1 F1 | P0/C1/R0 F1 | P1/C0/R1 F1 | R F1 | R PPV | R TPR | C F1 | C PPV | C TPR | P F1 | P PPV | P TPR |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| all | 18 | 68.918 | 78.957 | 68.571 | 78.637 | +10.039 | +30.521 | -6.533 | -0.347 | -3.523 | +4.375 | -0.320 | -0.010 | -0.463 |
| lower_duplication_fraction | 9 | 69.876 | 84.167 | 69.409 | 84.167 | +14.291 | +32.615 | -1.598 | -0.467 | -2.311 | +2.207 | +0.000 | +0.000 | +0.000 |
| missing_duplication_fraction | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| upper_duplication_fraction | 9 | 67.957 | 73.577 | 67.684 | 72.897 | +5.620 | +28.427 | -11.468 | -0.273 | -4.736 | +6.543 | -0.680 | -0.020 | -0.926 |

## Model Distance

| Bin | Families | P0/C0/R0 F1 | P0/C0/R1 F1 | P0/C1/R0 F1 | P1/C0/R1 F1 | R F1 | R PPV | R TPR | C F1 | C PPV | C TPR | P F1 | P PPV | P TPR |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| all | 18 | 68.918 | 78.957 | 68.571 | 78.637 | +10.039 | +30.521 | -6.533 | -0.347 | -3.523 | +4.375 | -0.320 | -0.010 | -0.463 |
| higher_than_median | 9 | 67.771 | 76.199 | 68.095 | 76.359 | +8.427 | +28.625 | -5.117 | +0.324 | -2.500 | +3.517 | +0.160 | +0.509 | +0.000 |
| lower_or_equal_median | 9 | 69.563 | 81.508 | 68.318 | 80.714 | +11.945 | +32.417 | -7.950 | -1.245 | -4.547 | +5.233 | -0.794 | -0.529 | -0.926 |

- Retrospective development-exposed fixed bins, not independent confirmation or causal effects.
- Initial HMM search stays on; P is downstream refinement, and R also changes pair semantics.
- Three conditional contrasts only, not a complete factorial or candidate/reconciliation interaction.
- Bins overlap and all-family rows repeat; neither is an independent finding.
- No subgroup intervals, significance, bootstrap, tuning, cutoffs or default promotion.
- Length/fragment flags do not establish fragment truth or completeness.
- Pfam descriptors are incomplete architecture annotations; reference duplication fractions are not ancestral histories.
- Estimated WAG+G4 distances include model/alignment/sampling/paralogy effects, not known genealogy or time.
- Direct retained-record/source checks, not repeated raw scoring or transitive scientific admission.
- Recovered cell timing stays failed/ineligible; shared-host timing distortion is unknown and tool-dependent.
