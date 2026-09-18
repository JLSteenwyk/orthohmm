# Corrected-Release QfO Scores

| Method | Status | GO similarity | EC similarity | VGNC F1 | SwissTrees F1 | TreeFam-A F1 | FAS | Secondary mean | Submitted pairs | Retained pairs |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| OrthoHMM high sensitivity | not_admitted | pending | pending | pending | pending | pending | pending | pending | pending | pending |
| OrthoHMM phylogeny satellite_v2 | not_admitted | pending | pending | pending | pending | pending | pending | pending | pending | pending |
| OrthoFinder 3.1.5 full | not_admitted | pending | pending | pending | pending | pending | pending | pending | pending | pending |
| OrthoFinder 3.1.5 sequence-only checkpoint | not_admitted | pending | pending | pending | pending | pending | pending | pending | pending | pending |
| SonicParanoid 2.0.9 | admitted | 0.454385 | 0.872761 | 0.982794 | 0.798459 | 0.771956 | 0.736680 | 0.769506 | 15248739 | 15248739 |
| ProteinOrtho 6.3.6 | admitted | 0.486336 | 0.963168 | 0.954896 | 0.718111 | 0.643187 | 0.813595 | 0.763216 | 4695385 | 4695385 |
| FastOMA 0.3.5 final orthologous groups | not_admitted | pending | pending | pending | pending | pending | pending | pending | pending | pending |
| OrthoMCL 1.4 | not_admitted | pending | pending | pending | pending | pending | pending | pending | pending | pending |

Only admitted corrected-release reports are included; pending is not zero and does not describe scheduler state.

GO/EC similarity and FAS are not F1; the six-metric mean is a project-defined secondary summary.

No paired uncertainty or corrected-release ranking is established by this point-estimate table.

This exporter checks admitted report hashes and native score arithmetic, not the complete inference/scorer workflow again.
