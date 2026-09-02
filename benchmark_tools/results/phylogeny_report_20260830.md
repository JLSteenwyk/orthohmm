# OrthoHMM phylogeny evaluation

## Decision summary

The conservative internally inferred-tree mode improves fresh OrthoBench
F-score from 70.4 to 70.9 by raising precision from 78.9 to 80.4 without
changing recall (63.5). It remains 1.8 F-score points below full OrthoFinder
3.1.5. On the frozen QfO assessment, its six-metric mean falls from 0.682548
to 0.662824 and remains below OrthoFinder's 0.782071. The implementation should
remain opt-in, and the current reconciliation rule should not become the
default: the small OrthoBench gain does not generalize and does not justify the
runtime cost.

The supplied-tree replay reaches F=73.0 and precision=86.2, but it uses the
species tree inferred by OrthoFinder and cached OrthoHMM gene trees. It shows
that a better species-tree topology can improve reconciliation, but is not a
fair win over OrthoFinder and is not an end-to-end runtime comparison.

## OrthoBench comparison

| Method | Tree information | Run type | F | Precision | Recall | Exact | Wall time | Peak RSS |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| OrthoFinder 3.1.5 sequence-only MCL | None | Checkpoint from full run | 58.7 | 45.7 | 82.0 | 18 | about 1,831 s to MCL | <=4.57 GiB run-wide maximum |
| OrthoHMM high sensitivity | None | Fresh end to end | 70.4 | 78.9 | 63.5 | 13 | 2,185.5 s | 11.90 GiB process tree |
| OrthoHMM supplied phylogeny | OrthoFinder species tree | Cached phylogeny-stage replay | 73.0 | 86.2 | 63.3 | 15 | 65.4 s stage only | 0.92 GiB stage process tree |
| OrthoHMM inferred phylogeny | Internal, 200 single-copy families | Fresh end to end | 70.9 | 80.4 | 63.5 | 14 | 3,023.4 s | 11.70 GiB process tree |
| OrthoFinder 3.1.5 full | Internal | Fresh end to end | 72.7 | 66.1 | 80.9 | 19 | 3,417.7 s | 4.57 GiB GNU-time maximum |

OrthoFinder's sequence-only time is derived from logged start and MCL
timestamps in the full run. Its stage-specific RSS was not measured, so the
full-run maximum is only an upper bound. GNU `time` reports the largest single
process for OrthoFinder, whereas the OrthoHMM harness samples the complete
process tree; those memory values are not strictly equivalent.

## QfO comparison

| Method | VGNC F | SwissTrees F | TreeFam-A F | EC | GO | FAS | Six-metric mean |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OrthoFinder 3.1.5 sequence-only | 0.14029994 | 0.69950736 | 0.71101320 | 0.75339446 | 0.40760020 | 0.55794138 | 0.54495942 |
| OrthoHMM high sensitivity | 0.66744809 | 0.67385086 | 0.57586996 | 0.93104489 | 0.47228231 | 0.77479433 | 0.68254841 |
| OrthoHMM inferred phylogeny | 0.74536115 | 0.64514895 | 0.34793251 | 0.97428495 | 0.48504497 | 0.77917375 | 0.66282438 |
| OrthoFinder 3.1.5 full | 0.98822936 | 0.85916044 | 0.74288699 | 0.94198234 | 0.46816557 | 0.69200063 | 0.78207089 |

The inferred-tree QfO candidate was frozen after independent synthetic
holdouts and OrthoBench. Native reconciled pairs were submitted directly; root
HOGs were not expanded into pairwise cliques. QfO was not used for tuning.

Relative to non-phylogenetic OrthoHMM, native pairs improve VGNC F by 0.0779,
EC by 0.0432, GO by 0.0128, and FAS by 0.0044, but reduce SwissTrees F by
0.0287 and TreeFam-A F by 0.2279. This is excessive conservatism: VGNC,
SwissTrees, and TreeFam-A precision rises to 0.9999, 0.9169, and 0.9662, while
recall falls to 0.5941, 0.4976, and 0.2122. TreeFam-A recall is the dominant
regression. Compared with full OrthoFinder, OrthoHMM wins EC, GO, and FAS but
loses all three reference-phylogeny F-scores and trails the mean by 0.1192.

The QfO inference used 976,504 genes and reconciled 26,159 of 390,657 candidate
families into 391,874 root HOGs and 2,659,577 native pairs. Instrumented wall
time was 74,011.3 seconds (20.56 hours) across the fresh run and verified
phylogeny restart, with 16.68 GiB peak process-tree RSS. The restart itself took
4,308.1 seconds and 3.34 GiB; QfO scoring took 1,587 seconds.

## What changed

The selected method does not merge candidate families. It reconciles 9,264
ambiguous families, bypasses 53,621 unambiguous families, and splits 364 of
62,885 source families into 63,514 root HOGs while preserving all 251,378
genes. Six split source families contain OrthoBench reference genes.

Post-inference analysis found that reconciliation added or removed no true
RefOG pair and introduced no cross-RefOG pair. Instead, it removed 191
reference-to-nonreference relations. All 191 came from multi-copy RefOGs;
162 were in the higher-identity half of RefOGs and 29 in the lower-identity
half. The modest gain is therefore paralogue/contaminant splitting in
multi-copy, relatively close families, not additional recovery of divergent
orthologues.

Permissive HMM candidate experiments did not generalize. Binary HMM graphs and
reciprocal family merging recovered recall but created weak transitive bridges
and severe over-grouping. Even broad OrthoFinder MCL candidates followed by
OrthoHMM reconciliation reached only F=68.8 with the conservative rule.
Relaxed root-HOG rules reached F=72.0-72.1 on that development input but
over-split a recent lineage-specific duplication in an independent holdout, so
they were rejected.

## Generalization and correctness

Three independent synthetic seeds cover ancient and recent duplication,
differential loss, remote homologues, fragments, multidomain proteins,
lineage-specific expansion, and composition bias. Supplied and inferred modes
recover all 32 root groups per seed with F=1.0. Inferred species trees have
rooted symmetric difference zero on every seed; native-pair F-scores are 1.0,
1.0, and 0.99175.

Production code contains no benchmark taxa, RefOG identifiers, or QfO labels.
Phylogeny remains disabled by default, and the non-phylogenetic path is
unchanged. Expensive alignments and trees are limited to ambiguous families,
are deterministic under the CPU budget, and have immutable restart
checkpoints. Missing optional programs and malformed or taxonomically invalid
trees fail with explicit errors.

The fresh QfO run exposed a general species-tree inference failure: one of 78
taxa occurred in no marker meeting the 50% occupancy cutoff, although it was
present in 304 lower-occupancy multi-species single-copy families. The fixed
selector retains 15 high-occupancy core markers and adds one connected
31-taxon placement marker, producing a 19,851-column supermatrix. The policy is
taxon-agnostic, requires one copy per represented species, records placement
markers in the checkpoint and manifest, and rejects disconnected markers. The
failed run and its verified restart are both retained in provenance.

## Outputs and provenance

The mode emits backward-compatible flat orthogroups, root and hierarchical
HOGs, native pairwise orthologues, raw/rooted/reconciled gene trees,
duplication/speciation records, a reconciliation summary, and a provenance
manifest. The fresh inferred OrthoBench root-HOG SHA-256 is
`b1911d1b9c31b2c72a95bb357a052cd83007cb5f9c36bb0f97afe1005e1399f0`,
identical to the frozen development replay.

Detailed machine-readable evidence is in:

- `benchmark_tools/results/phylogeny_development_20260830.json`
- `benchmark_tools/results/phylogeny_holdout_20260830.json`
- `benchmark_tools/results/candidate_superfamily_holdout_20260830.json`
- `benchmark_tools/results/orthobench_phylogeny_inferred_fresh_20260830.json`
- `benchmark_tools/results/qfo_phylogeny_inferred_20260902.json`
