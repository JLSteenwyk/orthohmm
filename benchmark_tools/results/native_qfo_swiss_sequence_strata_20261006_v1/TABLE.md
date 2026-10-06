# Native SwissTrees Frozen Sequence Strata

Descriptive percentages; F1 is the harmonic mean of macro precision and recall. Differences are R1 minus R0 in percentage points.
Empty bins are NA, not zero. These are development-exposed overlapping bins, with no new intervals or significance claims.

## all (18 families)

| Native cell | F1 | Precision | Recall |
|---|---:|---:|---:|
| p0_c0_r0 | 68.918 | 64.394 | 74.127 |
| p0_c0_r1 | 78.957 | 94.915 | 67.593 |
| R1 minus R0 | +10.039 | +30.521 | -6.533 |

## higher_entropy (9 families)

| Native cell | F1 | Precision | Recall |
|---|---:|---:|---:|
| p0_c0_r0 | 67.386 | 56.755 | 82.917 |
| p0_c0_r1 | 87.360 | 94.573 | 81.168 |
| R1 minus R0 | +19.974 | +37.818 | -1.749 |

## lower_entropy (9 families)

| Native cell | F1 | Precision | Recall |
|---|---:|---:|---:|
| p0_c0_r0 | 68.521 | 72.033 | 65.336 |
| p0_c0_r1 | 68.941 | 95.257 | 54.018 |
| R1 minus R0 | +0.420 | +23.224 | -11.318 |

## missing_entropy (0 families)

| Native cell | F1 | Precision | Recall |
|---|---:|---:|---:|
| p0_c0_r0 | NA | NA | NA |
| p0_c0_r1 | NA | NA | NA |
| R1 minus R0 | NA | NA | NA |

## concentrated (0 families)

| Native cell | F1 | Precision | Recall |
|---|---:|---:|---:|
| p0_c0_r0 | NA | NA | NA |
| p0_c0_r1 | NA | NA | NA |
| R1 minus R0 | NA | NA | NA |

## explicit_fragment (0 families)

| Native cell | F1 | Precision | Recall |
|---|---:|---:|---:|
| p0_c0_r0 | NA | NA | NA |
| p0_c0_r1 | NA | NA | NA |
| R1 minus R0 | NA | NA | NA |

## missing (0 families)

| Native cell | F1 | Precision | Recall |
|---|---:|---:|---:|
| p0_c0_r0 | NA | NA | NA |
| p0_c0_r1 | NA | NA | NA |
| R1 minus R0 | NA | NA | NA |

## no_explicit_fragment (18 families)

| Native cell | F1 | Precision | Recall |
|---|---:|---:|---:|
| p0_c0_r0 | 68.918 | 64.394 | 74.127 |
| p0_c0_r1 | 78.957 | 94.915 | 67.593 |
| R1 minus R0 | +10.039 | +30.521 | -6.533 |

## not_concentrated (18 families)

| Native cell | F1 | Precision | Recall |
|---|---:|---:|---:|
| p0_c0_r0 | 68.918 | 64.394 | 74.127 |
| p0_c0_r1 | 78.957 | 94.915 | 67.593 |
| R1 minus R0 | +10.039 | +30.521 | -6.533 |

## not_short_relative (11 families)

| Native cell | F1 | Precision | Recall |
|---|---:|---:|---:|
| p0_c0_r0 | 69.008 | 63.375 | 75.741 |
| p0_c0_r1 | 80.389 | 93.295 | 70.619 |
| R1 minus R0 | +11.381 | +29.920 | -5.121 |

## short_relative (7 families)

| Native cell | F1 | Precision | Recall |
|---|---:|---:|---:|
| p0_c0_r0 | 68.679 | 65.995 | 71.590 |
| p0_c0_r1 | 76.410 | 97.461 | 62.838 |
| R1 minus R0 | +7.731 | +31.465 | -8.753 |

Initial HMM search remains on in both cells. R0 uses group-clique pairs; R1 uses resolved native pairs. Global entropy is not divergence, relative shortness is not fragmentation, and absent fragment text does not prove completeness. Failed R1 timing stays ineligible.
