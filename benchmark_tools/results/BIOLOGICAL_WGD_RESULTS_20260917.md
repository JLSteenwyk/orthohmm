# WGD Biological Application

Development-exposed application; not orthology F1 or independent generalization.
All 240 experimental pairs retained; 239 input-eligible, 231 shared-pillar eligible.

| Method | Separated /239 | Supported /231 | Mean homolog coverage (%) | Unsupported splits /231 | Foreign-pillar cases /231 |
| --- | ---: | ---: | ---: | ---: | ---: |
| OrthoHMM high sensitivity | 58 | 56 | 99.149 | 0 | 56 |
| OrthoHMM phylogeny | 238 | 193 | 82.338 | 37 | 12 |
| OrthoFinder full (root HOGs) | 236 | 227 | 98.413 | 1 | 15 |
| SonicParanoid | 229 | 223 | 98.773 | 0 | 13 |
| OrthoFinder MCL checkpoint (diagnostic) | 190 | 184 | 99.221 | 0 | 20 |

Coverage can be high for merged groups. Foreign-pillar counts exclude unmapped members; neither is orthology precision.

## Paired Differences

All contrasts use 231 pairs. 20,000 paired pillar replicates; seed 20260920; differences in percentage points.

| First minus second | Endpoint | Difference | Bonferroni12 interval |
| --- | --- | ---: | --- |
| OrthoHMM phylogeny minus OrthoHMM high sensitivity | separation_rate | 75.325 | [66.667, 83.117] |
| OrthoHMM phylogeny minus OrthoHMM high sensitivity | supported_separation_rate | 59.307 | [49.784, 68.398] |
| OrthoHMM phylogeny minus OrthoHMM high sensitivity | mean_non_scer_coverage | -16.811 | [-20.563, -13.328] |
| OrthoHMM phylogeny minus OrthoFinder full (root HOGs) | separation_rate | 0.866 | [0.000, 3.030] |
| OrthoHMM phylogeny minus OrthoFinder full (root HOGs) | supported_separation_rate | -14.719 | [-21.645, -8.225] |
| OrthoHMM phylogeny minus OrthoFinder full (root HOGs) | mean_non_scer_coverage | -16.075 | [-19.755, -12.496] |
| OrthoHMM phylogeny minus SonicParanoid | separation_rate | 3.030 | [0.433, 6.926] |
| OrthoHMM phylogeny minus SonicParanoid | supported_separation_rate | -12.987 | [-19.913, -6.494] |
| OrthoHMM phylogeny minus SonicParanoid | mean_non_scer_coverage | -16.436 | [-20.101, -12.973] |
| OrthoHMM high sensitivity minus OrthoFinder full (root HOGs) | separation_rate | -74.459 | [-82.251, -65.801] |
| OrthoHMM high sensitivity minus OrthoFinder full (root HOGs) | supported_separation_rate | -74.026 | [-81.818, -64.935] |
| OrthoHMM high sensitivity minus OrthoFinder full (root HOGs) | mean_non_scer_coverage | 0.736 | [-0.524, 2.626] |

## Limitations

- Application uses development-exposed Saccharomyces data, not independent generalization.
- Homolog-supported paralog separation is not proof of cross-species copy-specific orthology.
- Foreign-pillar and unmapped members are separate descriptive diagnostics.
- MCL checkpoint is diagnostic, not a separately run sequence-only OrthoFinder pipeline.
- Caller must verify native admission and artifact identities before using this assembly.
