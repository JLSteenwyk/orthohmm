# External tool comparison (2026-09-04)

## QfO 2020

The QfO comparison uses the same retained six metrics for every method. The
reported mean is their unweighted arithmetic mean; it is a project summary,
not an official single QfO score.

| Rank | Method | VGNC F | SwissTrees F | TreeFam-A F | EC | GO | FAS | Mean |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | OrthoFinder 3.1.5 full | 0.988 | 0.859 | 0.743 | 0.942 | 0.468 | 0.692 | 0.782 |
| 2 | SonicParanoid 2.0.9 | 0.983 | 0.777 | 0.734 | 0.872 | 0.454 | 0.734 | 0.759 |
| 3 | ProteinOrtho 6.3.6 | 0.954 | 0.698 | 0.606 | 0.963 | 0.486 | 0.811 | 0.753 |
| 4 | OrthoHMM phylogeny `satellite_v2` | 0.901 | 0.792 | 0.575 | 0.969 | 0.490 | 0.762 | 0.748 |
| 5 | FastOMA 0.3.5 | 0.951 | 0.763 | 0.631 | 0.917 | 0.434 | 0.651 | 0.725 |
| 6 | OrthoHMM high sensitivity | 0.667 | 0.674 | 0.576 | 0.931 | 0.472 | 0.775 | 0.683 |
| 7 | OrthoFinder 3.1.5 sequence-only checkpoint | 0.140 | 0.700 | 0.711 | 0.753 | 0.408 | 0.558 | 0.545 |

SonicParanoid and ProteinOrtho are scored from their native pairwise outputs.
Their earlier group-clique means, 0.560 and 0.751, are invalid and retracted.

FastOMA produced 15,320,615 canonical pair rows, of which 15,277,489 passed
QfO identifier validation. Its retained-six mean is 0.724519. OrthoHMM
`satellite_v2` is higher by 0.023723, while full OrthoFinder is higher by
0.057551. FastOMA requires a supplied species tree; this run used the
OrthoFinder 3.1.5 tree inferred from the same 78 QfO proteomes. The tree has
exactly the expected 78 taxa and contains no benchmark labels.

## OrthoBench

The primary FastOMA result is its final `OrthologousGroups.tsv`. Its root-HOG
score is a useful diagnostic but is not substituted for the tool's final OG
output.

| Rank | Method | F-score (%) | Precision (%) | Recall (%) | Exact RefOGs | Native genes reported |
| ---: | --- | ---: | ---: | ---: | ---: | ---: |
| 1 | OrthoHMM phylogeny `satellite_v2` | 74.1 | 81.8 | 67.8 | 15 | 251,378 (100.0%) |
| 2 | OrthoFinder 3.1.5 full | 72.7 | 66.1 | 80.9 | 19 | 251,378 (100.0%) |
| 3 | OrthoHMM high sensitivity | 70.4 | 78.9 | 63.5 | 13 | 251,378 (100.0%) |
| 4 | OrthoFinder 3.1.5 sequence-only checkpoint | 58.7 | 45.7 | 82.0 | 18 | 251,378 (100.0%) |
| 5 | OrthoMCL 1.4 | 55.1 | 59.1 | 51.6 | 12 | 216,950 (86.3%) |
| 6 | SonicParanoid 2.0.9 | 46.8 | 36.5 | 64.9 | 19 | 206,466 (82.1%) |
| 7 | ProteinOrtho 6.3.6 | 45.1 | 97.1 | 29.3 | 6 | 184,897 (73.6%) |
| 8 | FastOMA 0.3.5 final OGs | 30.9 | 93.6 | 18.5 | 4 | 130,010 (51.7%) |

All eight rows were rerun through the official OrthoBench evaluator and also
recomputed with the repository's independent equal-RefOG implementation; both
calculations agree. “Native genes reported” counts genes in each tool's
retained native group output. The normalized SonicParanoid and ProteinOrtho
evaluator inputs add unassigned proteins as singleton groups to satisfy the
input-coverage requirement. Those singleton groups do not alter the pairwise
precision or recall.

On OrthoBench, OrthoHMM `satellite_v2` leads this set and exceeds full
OrthoFinder by 1.4 F-score points. On QfO, it ranks fourth and trails full
OrthoFinder by 0.033828 mean score. This split means there is no defensible
claim that one method wins “overall” without choosing an explicit weighting
between benchmarks.

Classic OrthoMCL QfO is running as inference job `20909`, followed by official
scoring job `20910`. It uses 180 BLAST threads and will submit validated native
cross-species ortholog/co-ortholog matrix edges. Both retained OrthoMCL
OrthoBench runs reproduce the reported score.

Exact values, source paths, job IDs, and checksums are in
`external_tool_comparison_20260904.json`.
