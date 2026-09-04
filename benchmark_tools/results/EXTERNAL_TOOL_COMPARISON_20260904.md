# External tool comparison (2026-09-04)

## QfO 2020

The QfO comparison uses the same retained six metrics for every method. The
reported mean is their unweighted arithmetic mean; it is a project summary,
not an official single QfO score.

| Rank | Method | VGNC F | SwissTrees F | TreeFam-A F | EC | GO | FAS | Mean |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | OrthoFinder 3.1.5 full | 0.988 | 0.859 | 0.743 | 0.942 | 0.468 | 0.692 | 0.782 |
| 2 | ProteinOrtho 6.3.6 | 0.893 | 0.743 | 0.642 | 0.948 | 0.480 | 0.800 | 0.751 |
| 3 | OrthoHMM phylogeny `satellite_v2` | 0.901 | 0.792 | 0.575 | 0.969 | 0.490 | 0.762 | 0.748 |
| 4 | FastOMA 0.3.5 | 0.951 | 0.763 | 0.631 | 0.917 | 0.434 | 0.651 | 0.725 |
| 5 | OrthoHMM high sensitivity | 0.667 | 0.674 | 0.576 | 0.931 | 0.472 | 0.775 | 0.683 |
| 6 | SonicParanoid 2.0.9 | 0.546 | 0.260 | 0.339 | 0.960 | 0.513 | 0.743 | 0.560 |
| 7 | OrthoFinder 3.1.5 sequence-only checkpoint | 0.140 | 0.700 | 0.711 | 0.753 | 0.408 | 0.558 | 0.545 |

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

| Rank | Method | F-score (%) | Precision (%) | Recall (%) | Exact RefOGs | Genes in OGs |
| ---: | --- | ---: | ---: | ---: | ---: | ---: |
| 1 | OrthoHMM phylogeny `satellite_v2` | 74.1 | 81.8 | 67.8 | 15 | 251,378 (100.0%) |
| 2 | OrthoFinder 3.1.5 full | 72.7 | 66.1 | 80.9 | 19 | 251,378 (100.0%) |
| 3 | OrthoHMM high sensitivity | 70.4 | 78.9 | 63.5 | 13 | 251,378 (100.0%) |
| 4 | OrthoFinder 3.1.5 sequence-only checkpoint | 58.7 | 45.7 | 82.0 | 18 | 251,378 (100.0%) |
| 5 | OrthoMCL 1.4 | 55.1 | 59.1 | 51.6 | 12 | 216,950 (86.3%) |
| 6 | SonicParanoid 2.0.9 | 46.8 | 36.5 | 64.9 | 19 | 206,466 (82.1%) |
| 7 | ProteinOrtho 6.3.6 | 45.1 | 97.1 | 29.3 | 6 | 184,897 (73.6%) |
| 8 | FastOMA 0.3.5 final OGs | 30.9 | 93.6 | 18.5 | 4 | 130,010 (51.7%) |

On OrthoBench, OrthoHMM `satellite_v2` leads this set and exceeds full
OrthoFinder by 1.4 F-score points. On QfO, it ranks third and trails full
OrthoFinder by 0.033828 mean score. This split means there is no defensible
claim that one method wins “overall” without choosing an explicit weighting
between benchmarks.

Classic OrthoMCL QfO was not launched. Its complete 12-proteome OrthoBench run
took 3 days 3 hours with a 9.7 GB all-versus-all BLAST table. Scaling that
legacy quadratic workflow to 78 QfO proteomes is not a practical comparator
run on the current system. Both retained OrthoMCL OrthoBench runs reproduce
the reported score.

Exact values, source paths, job IDs, and checksums are in
`external_tool_comparison_20260904.json`.
