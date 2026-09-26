# Corrected VGNC Count Decomposition

| Method | TP | FP | FN | Cross-block FP | Nonzero cross-block cells | F1 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| OrthoHMM high sensitivity | 19973 | 16078 | 3961 | 16076 | 7726 | 0.665933 |
| OrthoHMM phylogeny satellite_v2 | 19660 | 13 | 4274 | 13 | 8 | 0.901690 |
| OrthoFinder 3.1.5 full | 23519 | 130 | 415 | 130 | 38 | 0.988546 |
| OrthoFinder 3.1.5 sequence-only checkpoint | 23823 | 286004 | 111 | 286002 | 133852 | 0.142755 |
| SonicParanoid 2.0.9 | 23248 | 128 | 686 | 128 | 58 | 0.982794 |
| ProteinOrtho 6.3.6 | 22272 | 442 | 1662 | 442 | 168 | 0.954896 |
| FastOMA 0.3.5 final orthologous groups | 21742 | 92 | 2192 | 92 | 36 | 0.950096 |
| OrthoMCL 1.4 | 23520 | 25814 | 414 | 25814 | 12612 | 0.642027 |

Raw-row decomposition only; no new prediction-database rescore or biological confidence intervals.
OrthoFinder sequence-only is its diagnostic clustering checkpoint; FastOMA used a supplied tree.
