# Three Kingdoms: Audited Group Co-membership

| Method | Run | Precision | Recall | F1 | Reference coverage |
|---|---|---:|---:|---:|---:|
| OrthoHMM high sensitivity | historical | 1.000000 | 0.704026 | 0.826309 | 1.000000 |
| OrthoHMM phylogeny satellite_v2 | historical | 1.000000 | 0.773259 | 0.872133 | 1.000000 |
| OrthoFinder 3.1.5 full | historical | 1.000000 | 0.977421 | 0.988582 | 1.000000 |
| OrthoFinder 3.1.5 sequence-only checkpoint | historical | 0.990192 | 0.988711 | 0.989451 | 1.000000 |
| SonicParanoid 2.0.9 | contemporary matched input | 0.993443 | 0.989119 | 0.991276 | 0.998034 |
| ProteinOrtho 6.3.6 | historical | 1.000000 | 0.869967 | 0.930463 | 0.959705 |
| FastOMA 0.3.5 final orthologous groups | historical | 1.000000 | 0.756937 | 0.861655 | 0.894349 |
| OrthoMCL 1.4 | historical | 1.000000 | 0.978101 | 0.988929 | 0.996560 |
| SonicParanoid 2.0.9 (historical) | historical | 0.993436 | 0.988166 | 0.990794 | 0.998034 |

Historical Sonic is diagnostic only because of its known input mismatch.

- Supplementary conserved-family endpoint, not genome-wide orthology accuracy.
- Predictions outside the reference universe are not penalized.
- Historical normalized-group count audits do not prove identical historical input consumption.
- Historical Sonic is diagnostic only; the new run does not isolate the cause of its difference.
- No confidence intervals or superiority tests are supplied by this descriptive table.
- OrthoFinder sequence-only is an MCL checkpoint; FastOMA uses a supplied OrthoFinder tree.

## Input Evidence

- OrthoHMM high sensitivity: not_retained; historical consumption not proven.
- OrthoHMM phylogeny satellite_v2: all_12_staged_hashes_match; historical consumption not proven.
- OrthoFinder 3.1.5 full: all_12_staged_hashes_match; historical consumption not proven.
- OrthoFinder 3.1.5 sequence-only checkpoint: all_12_staged_hashes_match; historical consumption not proven.
- SonicParanoid 2.0.9: matched input verified.
- ProteinOrtho 6.3.6: not_retained; historical consumption not proven.
- FastOMA 0.3.5 final orthologous groups: all_12_staged_hashes_match; historical consumption not proven.
- OrthoMCL 1.4: all_12_staged_hashes_match; historical consumption not proven.
- SonicParanoid 2.0.9 (historical): known Danio input mismatch.
