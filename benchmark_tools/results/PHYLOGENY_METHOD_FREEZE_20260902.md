# Phylogeny method freeze, 2026-09-02

This record freezes the OrthoHMM phylogeny configuration before reading the
35-family OrthoBench validation partition or running a new QfO assessment.

## Frozen method

- High-sensitivity HMM seed groups with the bounded `satellite_v1` candidate
  expansion.
- Internally inferred species tree rooted by minimum root-to-tip variance over
  internal edges.
- Species-overlap root-HOG extraction.
- Positive-paralogy pair extraction, which removes pairs only at nodes with
  direct species-overlap evidence of duplication.
- One candidate expansion iteration with maximum satellite size 8, minimum
  source margin 1.5, minimum coverage 0.5, and minimum normalized support 0.02.

## Development evidence

| Method | F | Precision | Recall | Exact RefOGs |
|---|---:|---:|---:|---:|
| OrthoHMM high-sensitivity | 64.754 | 69.501 | 60.613 | 8 |
| Previous inferred phylogeny | 65.846 | 72.067 | 60.613 | 9 |
| OrthoFinder 3.1.5 full | 65.970 | 54.645 | 83.215 | 12 |
| Frozen OrthoHMM method | **68.197** | **71.216** | 65.424 | 9 |

The HMM candidate stage recovered 47.50% of development reference pairs
micro-averaged and 77.87% macro-averaged, with no candidate containing genes
from two known RefOGs. The selected method improves development F by 2.23
points over full OrthoFinder, while remaining substantially less sensitive.

## Negative results

- Midpoint rooting and unconstrained minimum-variance rooting both placed the
  root on the long *C. elegans* terminal branch.
- IQ-TREE 3.0.1 with `NQ.pfam+F+R4`, both with topology search and with a fixed
  topology, selected the same terminal root. The full search cost 1,099.6 s
  wall time and about 1.0 GB peak RSS.
- With midpoint rooting, the best satellite-candidate root rule reached only
  F=64.376 on development.

The exact configuration, checksums, ablations, and evidence paths are stored in
`phylogeny_method_freeze_20260902.json`.
