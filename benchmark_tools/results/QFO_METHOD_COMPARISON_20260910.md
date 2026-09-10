# Orthology Method Score Comparison (2026-09-10)

## QfO 2020

All completed methods below were evaluated by the same QfO 2020 pipeline. The
mean is this project's unweighted arithmetic mean of the six retained metrics;
it is not an official aggregate QfO score. Higher is better for every column.

| Rank | Method | VGNC F | SwissTrees F | TreeFam-A F | EC | GO | FAS | Mean | Status |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 1 | OrthoFinder 3.1.5 full | 0.9882 | 0.8592 | 0.7429 | 0.9420 | 0.4682 | 0.6920 | 0.7821 | Complete |
| 2 | SonicParanoid 2.0.9 | 0.9829 | 0.7771 | 0.7335 | 0.8717 | 0.4544 | 0.7338 | 0.7589 | Complete; native pairs |
| 3 | ProteinOrtho 6.3.6 | 0.9545 | 0.6978 | 0.6061 | 0.9628 | 0.4860 | 0.8113 | 0.7531 | Complete; native pairs |
| 4 | OrthoHMM phylogeny `satellite_v2` | 0.9007 | 0.7920 | 0.5755 | 0.9692 | 0.4901 | 0.7619 | 0.7482 | Complete |
| 5 | FastOMA 0.3.5 | 0.9512 | 0.7627 | 0.6307 | 0.9167 | 0.4345 | 0.6513 | 0.7245 | Complete |
| 6 | OrthoHMM high sensitivity | 0.6674 | 0.6739 | 0.5759 | 0.9310 | 0.4723 | 0.7748 | 0.6825 | Complete |
| 7 | OrthoFinder 3.1.5 sequence-only checkpoint | 0.1403 | 0.6995 | 0.7110 | 0.7534 | 0.4076 | 0.5579 | 0.5450 | Complete |
| - | OrthoMCL 1.4 | - | - | - | - | - | - | - | Running (`20909`; scoring `20910`) |

## Cross-Dataset Summary

This view uses the retained-six mean for QfO and F-score for OrthoBench and
Three Kingdoms. Scores are normalized to 0-1, but the benchmarks measure
different properties and should not be averaged without a declared weighting.

| Method | QfO mean | OrthoBench F | Three Kingdoms F |
| --- | ---: | ---: | ---: |
| OrthoFinder 3.1.5 full | 0.7821 | 0.7270 | 0.9886 |
| SonicParanoid 2.0.9 | 0.7589 | 0.4680 | 0.9908 |
| ProteinOrtho 6.3.6 | 0.7531 | 0.4510 | 0.9305 |
| OrthoHMM phylogeny `satellite_v2` | 0.7482 | 0.7411 | 0.8721 |
| FastOMA 0.3.5 | 0.7245 | 0.3090 | 0.8617 |
| OrthoHMM high sensitivity | 0.6825 | 0.7040 | 0.8263 |
| OrthoFinder 3.1.5 sequence-only checkpoint | 0.5450 | 0.5870 | 0.9895 |
| OrthoMCL 1.4 | Running | 0.5510 | 0.9889 |

The Three Kingdoms full-pipeline score restores OrthoFinder-sanitized
accessions from the exact staged FASTA headers before evaluation.

The corrected SonicParanoid and ProteinOrtho rows use native pairwise output;
their previous group-clique means of 0.5601 and 0.7511 are retracted.

The OrthoMCL QfO inference uses all 78 proteomes and 976,504 proteins with
native OrthoMCL 1.4 defaults, legacy NCBI BLAST 2.2.13, and MCL 02-063. Its
BLAST stage has 180 threads; the exact-compatible downstream pair stage has
64 checkpointed workers. The table will be updated from the official scoring
artifacts after jobs `20909` and `20910` complete. QfO pairs will be taken
from validated cross-species edges in OrthoMCL's native ortholog/co-ortholog
matrix, not from a clique expansion of its final clusters.
