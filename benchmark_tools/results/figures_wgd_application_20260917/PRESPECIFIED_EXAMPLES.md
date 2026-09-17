# All Prespecified WGD Examples

The six examples were selected by hash before native inference. None was replaced after seeing outcomes.
Support counts describe homologs in each anchor group, not validated ancestral-copy assignments.
High/Low/Sparse are source experimental classes, not orthology confidence labels.

## YER132C / YGL197W (High)

| Method | State | Support per anchor | Homolog coverage | Pillar groups | Foreign-pillar members |
| --- | --- | --- | --- | ---: | ---: |
| OrthoHMM high sensitivity | merged | 6, 6 | 6/6 | 1 | 0 |
| OrthoHMM phylogeny | separated | 3, 2 | 5/6 | 3 | 0 |
| OrthoFinder full | separated | 3, 3 | 6/6 | 2 | 0 |
| SonicParanoid | separated | 3, 3 | 6/6 | 2 | 0 |
| OrthoFinder MCL checkpoint (diagnostic) | separated | 3, 3 | 6/6 | 2 | 0 |

## YDR122W / YLR096W (High)

| Method | State | Support per anchor | Homolog coverage | Pillar groups | Foreign-pillar members |
| --- | --- | --- | --- | ---: | ---: |
| OrthoHMM high sensitivity | merged | 6, 6 | 6/6 | 1 | 0 |
| OrthoHMM phylogeny | separated | 3, 3 | 6/6 | 2 | 0 |
| OrthoFinder full | separated | 3, 3 | 6/6 | 2 | 0 |
| SonicParanoid | separated | 3, 3 | 6/6 | 2 | 0 |
| OrthoFinder MCL checkpoint (diagnostic) | merged | 6, 6 | 6/6 | 1 | 0 |

## YER059W / YIL050W (Low)

| Method | State | Support per anchor | Homolog coverage | Pillar groups | Foreign-pillar members |
| --- | --- | --- | --- | ---: | ---: |
| OrthoHMM high sensitivity | separated | 3, 3 | 6/6 | 2 | 0 |
| OrthoHMM phylogeny | separated | 3, 3 | 6/6 | 2 | 0 |
| OrthoFinder full | separated | 3, 3 | 6/6 | 2 | 0 |
| SonicParanoid | separated | 3, 3 | 6/6 | 2 | 0 |
| OrthoFinder MCL checkpoint (diagnostic) | separated | 3, 3 | 6/6 | 2 | 0 |

## YLR284C / YOR180C (Low)

| Method | State | Support per anchor | Homolog coverage | Pillar groups | Foreign-pillar members |
| --- | --- | --- | --- | ---: | ---: |
| OrthoHMM high sensitivity | merged | not evaluated | not evaluated | not evaluated | not evaluated |
| OrthoHMM phylogeny | separated | not evaluated | not evaluated | not evaluated | not evaluated |
| OrthoFinder full | separated | not evaluated | not evaluated | not evaluated | not evaluated |
| SonicParanoid | separated | not evaluated | not evaluated | not evaluated | not evaluated |
| OrthoFinder MCL checkpoint (diagnostic) | separated | not evaluated | not evaluated | not evaluated | not evaluated |

Reference exclusion retained: [{"lines": [7877, 1095], "reason": "anchors_in_different_pillars"}].

## YBR147W / YOL092W (Sparse)

| Method | State | Support per anchor | Homolog coverage | Pillar groups | Foreign-pillar members |
| --- | --- | --- | --- | ---: | ---: |
| OrthoHMM high sensitivity | merged | 6, 6 | 6/6 | 1 | 0 |
| OrthoHMM phylogeny | separated | 2, 3 | 5/6 | 3 | 0 |
| OrthoFinder full | separated | 3, 3 | 6/6 | 2 | 0 |
| SonicParanoid | separated | 3, 3 | 6/6 | 2 | 0 |
| OrthoFinder MCL checkpoint (diagnostic) | separated | 3, 3 | 6/6 | 2 | 0 |

## YCL048W / YDR522C (Sparse)

| Method | State | Support per anchor | Homolog coverage | Pillar groups | Foreign-pillar members |
| --- | --- | --- | --- | ---: | ---: |
| OrthoHMM high sensitivity | merged | 6, 6 | 6/6 | 1 | 7 |
| OrthoHMM phylogeny | separated | 0, 3 | 3/6 | 3 | 0 |
| OrthoFinder full | separated | 3, 3 | 6/6 | 2 | 0 |
| SonicParanoid | separated | 3, 3 | 6/6 | 2 | 0 |
| OrthoFinder MCL checkpoint (diagnostic) | separated | 3, 3 | 6/6 | 2 | 0 |

## Interpretation Limits

These tables establish final membership differences, not which search, tree or reconciliation decision caused them.
A pillar distributed over additional groups can reflect fragmentation, but the reference does not resolve correct ancestral copies.
See the complete machine-readable report for group IDs, sizes, unmapped members, missing assignments and all 240 pairs.
