# Three Kingdoms Retained-Source Audit

## Verified Scope

`audit_three_kingdoms_sources.py` traces the current retained compressed
downloads, raw FASTAs, staged FASTAs, BUSCO full tables and scored reference.
It does not download replacements or modify any input. Its source URLs are
read with a shell lexer from the checksum-pinned historical download script;
the script is never executed. Those URLs document download intent, not an
independent verification against original provider checksums.

All 12 compressed files decompress byte-for-byte to their retained raw
FASTAs. All input IDs are preserved, with no cross-species ID collisions.
Eleven staged FASTAs are byte-identical to raw copies. The zebrafish FASTA
has seven changed sequences, exactly explained by removing29 `*` markers.
The report retains the seven IDs, lengths and per-sequence stop counts.
None is in the scored reference. No other sequence-content change was found.

This is not proof of zero downstream effect, nor proof that every historical
tool consumed these exact staged bytes. Per-method input provenance remains
a separate check. No score was recomputed or changed by this audit.

## Lineage And Reference

All 12 retained BUSCO table headers report version5.8.2 and lineage
`eukaryota_odb10`, creation date2024-01-08,70 reference genomes and255 BUSCOs.
The retained `dataset.cfg` agrees and identifies OrthoDB10.1. Only Complete
single-copy hits enter the reference; Duplicated, Fragmented and Missing
entries are excluded. Groups require support in at least two species.

Reconstructing the reference independently from these tables reproduces the
scored reference byte-for-byte, SHA-256
`a5f3447056ecfa305442caff0d898d524eed350f13587d3247d3e1ba4c19757d`:
255groups,2,035genes,7,352within-group unordered pairs. Every retained Complete
hit maps to the corresponding species' staged FASTA. This validates the
reference construction from retained tables, not the biological correctness
of every BUSCO call or the original execution environment.

| Code | Staged proteins | Complete BUSCO families |
| --- | ---: | ---: |
| Amph | 43,615 | 226 |
| Anol | 19,176 | 221 |
| Arab | 48,321 | 108 |
| Asfu | 9,623 | 244 |
| Dani | 52,089 | 153 |
| Oryz | 42,582 | 162 |
| Scer | 6,600 | 237 |
| Scom | 13,194 | 240 |
| Sela | 34,825 | 27 |
| Spom | 5,146 | 243 |
| Xeno | 36,461 | 109 |
| Zeam | 131,585 | 65 |
| Total | 443,217 | Not a count of distinct families |

These are Complete-call counts, not a new accuracy endpoint or comparable
completeness fraction across proteomes with different isoform content.
Three Kingdoms remains a supplementary conserved-family benchmark.

## Acquisition And Rights

The pinned script specifies Ensembl release100 for Anolis/Danio,
Ensembl Metazoa47 for Amphimedon, Ensembl Plants47 for Arabidopsis/Selaginella/
Zea and58 for Oryza, and Ensembl Fungi47 for the four fungal proteomes.
Xenopus uses UniProt UP000186698 and a moving `current_release` URL. The
retained compressed checksum is therefore essential; rerunning that URL is
not a guarantee of recovering these input bytes. This audit does not assign
an unverified UniProt release date.

The current Ensembl disclaimer URL returned a service-unavailable page,
not usable licensing evidence. The [official September2025 archive](https://sep2025.archive.ensembl.org/info/about/legal/disclaimer.html)
states that project-generated data are unrestricted while warning of
third-party constraints. It distinguishes the code's Apache2.0 license from
data terms. This later archived policy does not itself clear every retained
historical proteome. Its snapshot is recorded in the data-rights register.
BUSCO lineage identity is now established; release-specific notices and the
scope of redistributed derived reference material still require review.

## Reproduce

```bash
python -m pytest -q tests/unit/test_audit_three_kingdoms_sources.py
python benchmark_tools/audit_three_kingdoms_sources.py --root three_kingdoms --output /tmp/three-kingdoms-source-audit.json
```

Use a new output path and the retained input tree.15 unit tests pass;
the full audit completed successfully with final hash rechecks. Its report
is `three_kingdoms_sources_20260918.json`. No inference or BUSCO rerun occurred.
