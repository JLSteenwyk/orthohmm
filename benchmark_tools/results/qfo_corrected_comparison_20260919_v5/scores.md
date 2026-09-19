# Corrected-Release QfO Publication Comparison

| Method | Status | GO similarity | EC similarity | VGNC F1 | SwissTrees F1 | TreeFam-A F1 | FAS | Secondary mean | Submitted pairs | Retained pairs | Mapping losses | Prediction semantics |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| OrthoHMM high sensitivity | admitted | 0.472271 | 0.932038 | 0.665933 | 0.685498 | 0.605008 | 0.776585 | 0.689555 | 9032719 | 9032719 | 0 | cross-species group-derived clique pairs |
| OrthoHMM phylogeny satellite_v2 | admitted | 0.490349 | 0.965650 | 0.901690 | 0.833513 | 0.614864 | 0.762993 | 0.761510 | 5959560 | 5959560 | 0 | native phylogenetically inferred pairs |
| OrthoFinder 3.1.5 full | admitted | 0.469548 | 0.936130 | 0.988546 | 0.848413 | 0.791918 | 0.691422 | 0.787663 | 14215382 | 14215382 | 0 | native phylogenetically inferred pairs |
| OrthoFinder 3.1.5 sequence-only checkpoint | admitted | 0.407043 | 0.753642 | 0.142755 | 0.690515 | 0.734407 | 0.561753 | 0.548352 | 163277439 | 163277439 | 0 | cross-species pre-phylogenetic MCL group-derived clique pairs (diagnostic) |
| SonicParanoid 2.0.9 | admitted | 0.454385 | 0.872761 | 0.982794 | 0.798459 | 0.771956 | 0.736680 | 0.769506 | 15248739 | 15248739 | 0 | native species-pair relations |
| ProteinOrtho 6.3.6 | admitted | 0.486336 | 0.963168 | 0.954896 | 0.718111 | 0.643187 | 0.813595 | 0.763216 | 4695385 | 4695385 | 0 | native post-clustering graph |
| FastOMA 0.3.5 final orthologous groups | not_admitted | not admitted | not admitted | not admitted | not admitted | not admitted | not admitted | not admitted | not admitted | not admitted | not admitted | native pairs; supplied OrthoFinder species tree |
| OrthoMCL 1.4 | not_admitted | not admitted | not admitted | not admitted | not admitted | not admitted | not admitted | not admitted | not admitted | not admitted | not admitted | final MCL group-derived pairs |

Only supplied corrected-release admissions are included; missing is not zero or a scheduler-state claim.

OrthoHMM high sensitivity is p1_c0_r0; phylogeny satellite_v2 is p1_c1_r1. Other factorial cells are ablations, not interchangeable defaults.

The pinned independent replay admission establishes equality with the corrected native high-sensitivity partition, not historical-input outputs.

GO/EC similarity and FAS are not F1; the six-metric mean is a project-defined secondary summary.

Prediction semantics differ by method and are shown explicitly. Pair volume is not protein coverage.

This table establishes no paired significance, ranking, independent generalization or matched efficiency.

Checks cover report hashes, conversion binding and native score arithmetic; upstream admission audits are not rerun.
