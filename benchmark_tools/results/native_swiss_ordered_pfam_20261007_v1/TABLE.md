# Ordered-Pfam Annotation Strata

Descriptive development-exposed points; no new intervals.

| Stratum | Families | Cell | F1 (%) | Precision (%) | Recall (%) |
| --- | ---: | --- | ---: | ---: | ---: |
| all | 18 | p0_c0_r0 | 68.918 | 64.394 | 74.127 |
| all | 18 | p0_c0_r1 | 78.957 | 94.915 | 67.593 |
| all | 18 | p0_c1_r0 | 68.571 | 60.871 | 78.502 |
| all_members_usable_same_signature | 4 | p0_c0_r0 | 74.632 | 66.532 | 84.978 |
| all_members_usable_same_signature | 4 | p0_c0_r1 | 81.645 | 95.997 | 71.027 |
| all_members_usable_same_signature | 4 | p0_c1_r0 | 73.066 | 61.578 | 89.824 |
| all_members_usable_multiple_signatures | 1 | p0_c0_r0 | 65.327 | 49.242 | 97.015 |
| all_members_usable_multiple_signatures | 1 | p0_c0_r1 | 97.015 | 97.015 | 97.015 |
| all_members_usable_multiple_signatures | 1 | p0_c1_r0 | 65.327 | 49.242 | 97.015 |
| some_members_unusable | 13 | p0_c0_r0 | 66.901 | 64.902 | 69.027 |
| some_members_unusable | 13 | p0_c0_r1 | 76.484 | 94.421 | 64.274 |
| some_members_unusable | 13 | p0_c1_r0 | 67.034 | 61.547 | 73.594 |

| Stratum | Families | Contrast | Delta F1 (pp) | Delta Precision (pp) | Delta Recall (pp) |
| --- | ---: | --- | ---: | ---: | ---: |
| all | 18 | R_at_P0_C0 | +10.039 | +30.521 | -6.533 |
| all | 18 | C_at_P0_R0 | -0.347 | -3.523 | +4.375 |
| all_members_usable_same_signature | 4 | R_at_P0_C0 | +7.013 | +29.464 | -13.951 |
| all_members_usable_same_signature | 4 | C_at_P0_R0 | -1.566 | -4.954 | +4.846 |
| all_members_usable_multiple_signatures | 1 | R_at_P0_C0 | +31.688 | +47.773 | +0.000 |
| all_members_usable_multiple_signatures | 1 | C_at_P0_R0 | +0.000 | +0.000 | +0.000 |
| some_members_unusable | 13 | R_at_P0_C0 | +9.583 | +29.519 | -4.754 |
| some_members_unusable | 13 | C_at_P0_R0 | +0.133 | -3.354 | +4.567 |

Initial HMM search on, downstream profiles off; conditional C/R effects, not an interaction.
Predicted domain order and unusable annotations do not establish complete architecture or biological absence.
Failed R1 timing stays ineligible; no raw scoring, bootstrap, inference, default or independent-validation claim.
