# All Methods Across All Benchmark Datasets (2026-09-10)

All values use a 0-1 scale and higher is better. QfO is the unweighted mean
of the six retained QfO metrics (VGNC F, SwissTrees F, TreeFam-A F, EC, GO,
and FAS). OrthoBench and Three Kingdoms report F1 scores from their respective
common evaluators.

| Method | QfO mean | OrthoBench F1 | Three Kingdoms F1 | Complete datasets |
| --- | ---: | ---: | ---: | ---: |
| OrthoFinder 3.1.5 full phylogenetic pipeline | **0.7821** | 0.7270 | 0.8561 | 3/3 |
| ProteinOrtho 6.3.6 | 0.7511 | 0.4510 | 0.9305 | 3/3 |
| OrthoHMM 0.5.0 phylogeny `satellite_v2` | 0.7482 | **0.7411** | 0.8721 | 3/3 |
| FastOMA 0.3.5 | 0.7245 | 0.3090 | 0.8617 | 3/3 |
| OrthoHMM 0.5.0 high sensitivity | 0.6825 | 0.7040 | 0.8263 | 3/3 |
| SonicParanoid 2.0.9 | 0.5601 | 0.4680 | 0.8859 | 3/3 |
| OrthoFinder 3.1.5 sequence-only checkpoint | 0.5450 | 0.5870 | **0.9895** | 3/3 |
| OrthoMCL 1.4 | Running | 0.5510 | 0.9889 | 2/3 |

The best completed QfO value is OrthoFinder full at **0.7821**. The best
OrthoBench F1 is OrthoHMM phylogeny at **0.7411**, and the best Three Kingdoms
F1 is OrthoFinder sequence-only at **0.9895**. OrthoMCL QfO inference job
`20904` is running, with official QfO scoring job `20905` queued after it.

No cross-dataset average is reported because the QfO mean and the two F1
scores measure different benchmark objectives. An aggregate would require an
explicit, scientifically justified weighting across datasets.
