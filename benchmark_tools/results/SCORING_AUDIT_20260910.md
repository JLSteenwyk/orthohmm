# Benchmark Scoring Audit (2026-09-10)

This audit covers the retained cross-tool comparisons on QfO 2020,
OrthoBench, and Three Kingdoms, plus historical external scores referenced by
the OrthoHMM hill-climb reports. It validates score arithmetic, identifier
mapping, group partitions, pair formatting, and whether native pairwise output
was used when a tool provides it.

## Outcome

| Benchmark | Result |
| --- | --- |
| OrthoBench | All eight published rows reproduce with the official evaluator and an independent equal-RefOG implementation. |
| Three Kingdoms | All valid rows reproduce from TP/FP/FN counts. SonicParanoid and historical OrthoFinder 2.5.5 rows required corrections. |
| QfO 2020 | The six metric calculations and unweighted project mean reproduce from official JSON. SonicParanoid and ProteinOrtho required corrected native-pair runs; OrthoMCL remains in progress. |

## Corrections

| Method and benchmark | Retracted value | Correct value | Cause |
| --- | ---: | ---: | --- |
| SonicParanoid 2.0.9, Three Kingdoms F1 | 0.8859 | 0.990794 | Metadata columns were treated as genes and comma-delimited species cells were not split. |
| OrthoFinder 2.5.5 DIAMOND, Three Kingdoms F1 | 0.8552 | 0.987885 | Sanitized identifiers were not restored from the staged FASTA headers. |
| OrthoFinder 2.5.5 MMseqs, Three Kingdoms F1 | 0.8518 | 0.983899 | Sanitized identifiers were not restored from the staged FASTA headers. |
| SonicParanoid 2.0.9, QfO mean | 0.560060 | 0.758914 | QfO pairs were clique-expanded from the same malformed global-group conversion instead of using native species-pair relations. |
| ProteinOrtho 6.3.6, QfO mean | 0.751068 | 0.753099 | The old score clique-expanded global groups even though ProteinOrtho emits a native pairwise orthology graph. |

The corrections do not reverse an OrthoHMM acceptance decision. They do
invalidate the earlier external-tool rankings and strengthen the conclusion
that the historical OrthoHMM high-sensitivity configuration did not lead the
Three Kingdoms benchmark.

## OrthoBench

| Method | F-score (%) | Precision (%) | Recall (%) | Exact RefOGs |
| --- | ---: | ---: | ---: | ---: |
| OrthoHMM phylogeny `satellite_v2` | 74.106074 | 81.770454 | 67.755336 | 15 |
| OrthoFinder 3.1.5 full | 72.736480 | 66.065103 | 80.906577 | 19 |
| OrthoHMM high sensitivity | 70.358998 | 78.949544 | 63.454479 | 13 |
| OrthoFinder 3.1.5 sequence-only | 58.705963 | 45.706163 | 82.039827 | 18 |
| OrthoMCL 1.4 | 55.065342 | 59.078055 | 51.563064 | 12 |
| SonicParanoid 2.0.9 | 46.757609 | 36.535180 | 64.922809 | 19 |
| ProteinOrtho 6.3.6 | 45.057347 | 97.056208 | 29.338789 | 6 |
| FastOMA 0.3.5 final OGs | 30.906942 | 93.631476 | 18.508164 | 4 |

The report's gene-count column now means genes present in the retained native
tool output. SonicParanoid and ProteinOrtho inputs were padded with unassigned
singletons solely to meet the official evaluator's coverage requirement;
singletons do not change pairwise precision or recall.

## Three Kingdoms

The scorer uses micro-averaged pairs among the 2,035 genes in 255 BUSCO
reference groups. Exact values below are recomputed from TP, FP, and FN rather
than parsed from rounded display text.

| Method | TP | FP | FN | F1 |
| --- | ---: | ---: | ---: | ---: |
| SonicParanoid 2.0.9 | 7,265 | 48 | 87 | 0.990794 |
| OrthoFinder 3.1.5 sequence-only | 7,269 | 72 | 83 | 0.989451 |
| OrthoMCL 1.4 | 7,191 | 0 | 161 | 0.988929 |
| OrthoFinder 3.1.5 full | 7,186 | 0 | 166 | 0.988582 |
| ProteinOrtho 6.3.6 | 6,396 | 0 | 956 | 0.930463 |
| OrthoHMM phylogeny `satellite_v2` | 5,685 | 0 | 1,667 | 0.872133 |
| FastOMA 0.3.5 final OGs | 5,565 | 0 | 1,787 | 0.861655 |
| OrthoHMM high sensitivity | 5,176 | 0 | 2,176 | 0.826309 |

This is a narrow secondary metric. It cannot count false-positive relations
involving genes outside the BUSCO reference set, which explains why several
methods have precision 1.0 despite materially different full-dataset outputs.

## QfO Validation

All retained pair files were checked for exactly two tab-separated columns,
canonical ordering, self-pairs, duplicate rows, and QfO identifier validity.
Group-derived methods were also checked for duplicate or unknown genes.
Native relations are used for OrthoFinder full, OrthoHMM phylogeny, FastOMA,
SonicParanoid, and ProteinOrtho. Clique expansion remains appropriate for the
group-only OrthoHMM high-sensitivity and OrthoFinder sequence-only outputs.
The pending OrthoMCL assessment will use its cross-species native
ortholog/co-ortholog matrix edges; expanding final MCL groups would introduce
cross-species pairs that OrthoMCL did not call directly.

| Method | VGNC F | SwissTrees F | TreeFam-A F | EC | GO | FAS | Mean |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OrthoFinder 3.1.5 full | 0.988229 | 0.859160 | 0.742887 | 0.941982 | 0.468166 | 0.692001 | 0.782071 |
| SonicParanoid 2.0.9, native pairs | 0.982942 | 0.777107 | 0.733531 | 0.871666 | 0.454413 | 0.733824 | 0.758914 |
| ProteinOrtho 6.3.6, native pairs | 0.954497 | 0.697835 | 0.606133 | 0.962838 | 0.485962 | 0.811332 | 0.753099 |
| OrthoHMM phylogeny `satellite_v2` | 0.900730 | 0.791976 | 0.575495 | 0.969195 | 0.490120 | 0.761939 | 0.748243 |
| FastOMA 0.3.5 | 0.951223 | 0.762728 | 0.630710 | 0.916692 | 0.434480 | 0.651283 | 0.724519 |
| OrthoHMM high sensitivity | 0.667448 | 0.673851 | 0.575870 | 0.931045 | 0.472282 | 0.774794 | 0.682548 |
| OrthoFinder 3.1.5 sequence-only | 0.140300 | 0.699507 | 0.711013 | 0.753394 | 0.407600 | 0.557941 | 0.544959 |

Corrected SonicParanoid provenance:

- 3,003 native species-pair tables, the complete matrix for 78 species
- 15,022,677 distinct native pairs after removing 5,715 duplicate relations
- 14,975,235 QfO-valid pairs after removing 47,442 unmapped relations
- native pairs SHA-256: `734b99a33036c47d9ef0a12e1331a317682c991fa5f26b2fdeafb0427bd95fb6`
- filtered pairs SHA-256: `ce01ac6492c3d59fe79e2b4298d258d864ea2925b93a02e93d257f927f4c0741`
- official scoring job `20907`, completed successfully in 00:30:52

Corrected ProteinOrtho provenance:

- 4,579,157 native relations from `qfo.proteinortho-graph`
- 4,563,756 QfO-valid pairs after removing 15,401 unmapped relations
- native pairs SHA-256: `7d72e0c6778544813e7a4d7518f1983447458dc592fbed99da903601c1ce8fae`
- filtered pairs SHA-256: `2245d37f7a35e8764dcd35adc914a21e5f5b31fb9791478e5fc5a78445285981`
- official scoring job `20908`, completed successfully in 00:24:07
- official result JSON SHA-256 values: VGNC `fbe592028d0a02fdbc9729517cabcdd2337177179804db3c4b3d57d5df88f529`, SwissTrees `eb26d9efeb37a3d79799522eaa6771160b960009d0d139956c22657a9e04ba52`, TreeFam-A `bd1d8504916f9abfe2ad6ace629a038f7586d32c6ae3848fe12bcae4aa706dcc`, EC `1019637c0bcc785b12bcd193d0633c88b71bcda86e78d6a39dce2551f1a2a3e7`, GO `8d0defa00e949b07f09e2a687cf6820573ccf52932586abffdd4d0d85b389dc9`, and FAS `8d857e991d222377ec0ae93b6223ff53ad0193ef78d5292cfb41b325af13c753`

The QfO “mean” is an explicitly project-defined, unweighted arithmetic mean
of VGNC F1, SwissTrees F1, TreeFam-A F1, EC, GO, and FAS. Its arithmetic is
correct, but it is not an official QfO aggregate and should not be interpreted
as one.

## Guardrails Added

- Strict SonicParanoid metadata, cell, count, and cross-group validation
- Native SonicParanoid and ProteinOrtho QfO pair converters
- Native OrthoMCL matrix-edge QfO pair converter with reciprocal-score checks
- Complete and unique species-pair matrix validation
- Canonical pair, self-pair, duplicate-gene, and unknown-gene rejection
- Exact Three Kingdoms metrics derived from TP/FP/FN with consistency checks
- Regression tests for each malformed-input case

OrthoMCL 1.4 QfO inference job `20909` and dependent scoring job `20910` are
still running and are not assigned a QfO score in this audit. Jobs `20904` and
`20905` were canceled before producing a result when this audit identified the
cluster-clique conversion; only 1.6% of the BLAST query set had been processed.
