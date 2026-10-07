# Native SwissTrees Fixed Strata: Three Cells

F1 in percent; changes in percentage points against p0_c0_r0. NA denotes an empty bin.
R: reconciliation at P0/C0. C: candidate expansion at P0/R0. Initial HMM search on in all cells.
These are development-exposed descriptions, not intervals or evidence of a causal interaction.

## Sequence

| Bin | Families | Baseline F1 | R1 F1 | C1 F1 | R F1 change | R precision change | R recall change | C F1 change | C precision change | C recall change |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| all | 18 | 68.918 | 78.957 | 68.571 | +10.039 | +30.521 | -6.533 | -0.347 | -3.523 | +4.375 |
| concentrated | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| explicit_fragment | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| higher_entropy | 9 | 67.386 | 87.360 | 66.650 | +19.974 | +37.818 | -1.749 | -0.736 | -2.345 | +3.079 |
| lower_entropy | 9 | 68.521 | 68.941 | 69.120 | +0.420 | +23.224 | -11.318 | +0.599 | -4.702 | +5.671 |
| missing | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| missing_entropy | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| no_explicit_fragment | 18 | 68.918 | 78.957 | 68.571 | +10.039 | +30.521 | -6.533 | -0.347 | -3.523 | +4.375 |
| not_concentrated | 18 | 68.918 | 78.957 | 68.571 | +10.039 | +30.521 | -6.533 | -0.347 | -3.523 | +4.375 |
| not_short_relative | 11 | 69.008 | 80.389 | 68.772 | +11.381 | +29.920 | -5.121 | -0.236 | -3.961 | +5.889 |
| short_relative | 7 | 68.679 | 76.410 | 67.976 | +7.731 | +31.465 | -8.753 | -0.703 | -2.836 | +1.996 |

## Domain

| Bin | Families | Baseline F1 | R1 F1 | C1 F1 | R F1 change | R precision change | R recall change | C F1 change | C precision change | C recall change |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| all | 18 | 68.918 | 78.957 | 68.571 | +10.039 | +30.521 | -6.533 | -0.347 | -3.523 | +4.375 |
| median_pfam_types_at_least_two | 6 | 63.986 | 84.932 | 64.629 | +20.946 | +42.987 | -1.262 | +0.643 | -0.958 | +3.839 |
| median_pfam_types_below_two | 12 | 71.006 | 75.897 | 70.334 | +4.891 | +24.288 | -9.169 | -0.672 | -4.806 | +4.643 |
| repeated_type_fraction_at_least_quarter | 3 | 65.595 | 74.731 | 66.582 | +9.136 | +37.390 | -13.296 | +0.987 | -1.054 | +4.496 |
| repeated_type_fraction_below_quarter | 15 | 69.545 | 79.749 | 68.954 | +10.204 | +29.147 | -5.181 | -0.591 | -4.017 | +4.351 |

## Duplication

| Bin | Families | Baseline F1 | R1 F1 | C1 F1 | R F1 change | R precision change | R recall change | C F1 change | C precision change | C recall change |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| all | 18 | 68.918 | 78.957 | 68.571 | +10.039 | +30.521 | -6.533 | -0.347 | -3.523 | +4.375 |
| lower_duplication_fraction | 9 | 69.876 | 84.167 | 69.409 | +14.291 | +32.615 | -1.598 | -0.467 | -2.311 | +2.207 |
| missing_duplication_fraction | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| upper_duplication_fraction | 9 | 67.957 | 73.577 | 67.684 | +5.620 | +28.427 | -11.468 | -0.273 | -4.736 | +6.543 |

- Development-exposed descriptive fixed bins; no subgroup intervals, significance or tuning.
- Macro family precision/recall then harmonic F1; not pooled pairs or mean-family F1.
- Initial HMM search on and downstream profile refinement off in all three cells.
- Two conditional contrasts against baseline; candidate/reconciliation interaction not identified.
- Bins overlap and all-family rows repeat; they are not independent findings.
- Length/composition proxies are not calibrated divergence or literal fragment truth.
- No explicit fragment annotation is not evidence of completeness.
- Pfam descriptors are not complete validated domain architectures or domain-loss truth.
- Reference informative-node duplication fractions are not ancestral duplication histories; default S is not explicit speciation.
- Direct artifact/source checks only; prior raw validation inherited, not repeated or transitively readmitted.
- Cell 7 failed timing remains ineligible; no timing repair or new inference resources.
- Shared-host contention effects are unknown and potentially tool-dependent, not isolated speed evidence.
