### Profile Refinement Across All Fixed Bins

All 23 conditional profile-bin differences are reported below in percentage
points. Overlapping and repeated all-family bins are audit rows, not independent
findings. NA denotes an empty bin, not zero. No subgroup intervals were computed.

| Suite | Bin | Families | F1 Difference | PPV Difference | TPR Difference |
| --- | --- | ---: | ---: | ---: | ---: |
| sequence | all | 18 | -0.320 | -0.010 | -0.463 |
| sequence | concentrated | 0 | NA | NA | NA |
| sequence | explicit_fragment | 0 | NA | NA | NA |
| sequence | higher_entropy | 9 | -0.763 | -0.529 | -0.926 |
| sequence | lower_entropy | 9 | +0.133 | +0.509 | +0.000 |
| sequence | missing | 0 | NA | NA | NA |
| sequence | missing_entropy | 0 | NA | NA | NA |
| sequence | no_explicit_fragment | 18 | -0.320 | -0.010 | -0.463 |
| sequence | not_concentrated | 18 | -0.320 | -0.010 | -0.463 |
| sequence | not_short_relative | 11 | -0.499 | -0.016 | -0.758 |
| sequence | short_relative | 7 | +0.000 | +0.000 | +0.000 |
| domain | all | 18 | -0.320 | -0.010 | -0.463 |
| domain | median_pfam_types_at_least_two | 6 | +0.000 | +0.000 | +0.000 |
| domain | median_pfam_types_below_two | 12 | -0.496 | -0.015 | -0.694 |
| domain | repeated_type_fraction_at_least_quarter | 3 | +0.000 | +0.000 | +0.000 |
| domain | repeated_type_fraction_below_quarter | 15 | -0.378 | -0.012 | -0.556 |
| duplication | all | 18 | -0.320 | -0.010 | -0.463 |
| duplication | lower_duplication_fraction | 9 | +0.000 | +0.000 | +0.000 |
| duplication | missing_duplication_fraction | 0 | NA | NA | NA |
| duplication | upper_duplication_fraction | 9 | -0.680 | -0.020 | -0.926 |
| model_distance | all | 18 | -0.320 | -0.010 | -0.463 |
| model_distance | higher_than_median | 9 | +0.160 | +0.509 | +0.000 |
| model_distance | lower_or_equal_median | 9 | -0.794 | -0.529 | -0.926 |

Overall F1 changes by -0.320215 percentage points.
Lower-entropy families have +0.132851 points, while higher-entropy
families have -0.762755. Above-median model-distance families
have +0.160465, versus -0.794278
at or below the fixed median. Each entropy/distance bin has nine families.

These patterns reflect the same retained changed families rather than many
independent successes. Profile-minus-reference integer count changes are:

| Family | TP | FP | FN | TN |
| --- | ---: | ---: | ---: | ---: |
| CASP | +0 | -3 | +0 | +3 |
| GH14 | -1 | +0 | +1 | +0 |

The other 16 families have unchanged counts. CASP contributes
the precision gain in the lower-entropy/higher-distance bins; GH14 contributes
the recall loss in the complementary bins. Both occur in the not-short-relative,
lower-Pfam-type, lower-repeat-fraction and upper-reference-duplication bins.
Their complementary bins are unchanged. Five empty bins support no effect.
The [existing complete pair localization](NATIVE_PROFILE_LOCALIZATION_RESULT_20261007.md)
places all four removed calls at pre-reconciliation candidate separation.
It does not identify a causal HMM edge or establish correct gene/species trees.

These are development-exposed descriptive associations, not a general entropy
or divergence benefit, an independent validation, or an initial-HMM effect.
Existing aggregate intervals are not subgroup intervals. Fragment flags, Pfam
descriptors, reference duplication fractions and model distances retain their
proxy limitations. Missing factorial cells and failed/ineligible timing remain.
[Complete four-cell table](native_qfo_four_cell_strata_20261009_v2/TABLE.md),
[full-precision result](native_qfo_four_cell_strata_20261009_v2/report.json),
[actual commands/readback](native_qfo_four_cell_strata_execution_20261009_v2.json).

