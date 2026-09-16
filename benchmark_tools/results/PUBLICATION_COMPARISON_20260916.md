# Retained Benchmark Comparison

Work in progress; not a frozen publication baseline.

| Method | OrthoBench F1 (%) | QfO six-metric mean | Three Kingdoms F1 |
| --- | ---: | ---: | ---: |
| OrthoHMM high sensitivity | 70.358998 | 0.682548 | 0.826309 |
| OrthoHMM phylogeny satellite_v2 | 74.106074 | 0.748243 | 0.872133 |
| OrthoFinder 3.1.5 full | 72.736480 | 0.782071 | 0.988582 |
| OrthoFinder 3.1.5 sequence-only checkpoint | 58.705963 | 0.544959 | 0.989451 |
| SonicParanoid 2.0.9 | 46.757609 | 0.758914 | 0.990794 |
| ProteinOrtho 6.3.6 | 45.057347 | 0.753099 | 0.930463 |
| FastOMA 0.3.5 final orthologous groups | 30.906942 | 0.724519 | 0.861655 |
| OrthoMCL 1.4 | 55.065342 | pending | 0.988929 |

## QfO Components

| Method | VGNC F | SwissTrees F | TreeFam-A F | EC | GO | FAS | Output |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| OrthoHMM high sensitivity | 0.667448 | 0.673851 | 0.575870 | 0.931045 | 0.472282 | 0.774794 | group-derived pairs |
| OrthoHMM phylogeny satellite_v2 | 0.900730 | 0.791976 | 0.575495 | 0.969195 | 0.490120 | 0.761939 | phylogenetically inferred pairs |
| OrthoFinder 3.1.5 full | 0.988229 | 0.859160 | 0.742887 | 0.941982 | 0.468166 | 0.692001 | phylogenetically inferred pairs |
| OrthoFinder 3.1.5 sequence-only checkpoint | 0.140300 | 0.699507 | 0.711013 | 0.753394 | 0.407600 | 0.557941 | MCL checkpoint group-derived pairs |
| SonicParanoid 2.0.9 | 0.982942 | 0.777107 | 0.733531 | 0.871666 | 0.454413 | 0.733824 | native species-pair relations |
| ProteinOrtho 6.3.6 | 0.954497 | 0.697835 | 0.606133 | 0.962838 | 0.485962 | 0.811332 | native post-clustering graph |
| FastOMA 0.3.5 final orthologous groups | 0.951223 | 0.762728 | 0.630710 | 0.916692 | 0.434480 | 0.651283 | native pairs; supplied OrthoFinder species tree |
| OrthoMCL 1.4 | pending | pending | pending | pending | pending | pending | final MCL group-derived pairs |

OrthoMCL pre-clustering diagnostic mean: 0.724414. This is not the final-group assessment and is excluded from the main comparison.

## Limitations

- QfO and OrthoBench are development-exposed; these results do not establish independent generalization.
- QfO mean is a project-defined unweighted six-metric secondary summary, not an official score.
- QfO metric_x means precision/recall or assessed relation count depending on the recorded axis, not total prediction coverage.
- Three Kingdoms scores only pairs among BUSCO reference genes, ignoring false positives involving other genes.
- OrthoBench non-primary-tool rows are from the prior scoring audit; their raw-output provenance is not yet consolidated here.
- Runtime budgets and memory accounting differ. Listed Three Kingdoms timings are historical evidence, not a matched efficiency experiment.
- FastOMA used supplied trees; tree-construction cost and reference-resource overlap require explicit accounting.
- Independent validation, controlled ablations, error strata, robustness, biological application, and release work remain open.
