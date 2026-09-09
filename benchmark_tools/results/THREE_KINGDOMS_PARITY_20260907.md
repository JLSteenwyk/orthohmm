# Three Kingdoms Method Parity Benchmark

All methods were evaluated against the same 255 BUSCO reference orthogroups (2,035 reference genes) from the same 12-proteome, 443,217-protein dataset.

| Method | Variant | Precision | Recall | F-score | Ref. coverage | Predicted OGs | Wall time | Peak RSS |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| OrthoFinder 3.1.5 | sequence-only MCL checkpoint | 0.9902 | 0.9887 | 0.9895 | 100.0% | 65,602 | 1:21:24 | n/a |
| OrthoMCL 1.4 | defaults; legacy BLASTP; byte-compatible conversion and pair parallelism; MCL inflation 1.5 | 1.0000 | 0.9781 | 0.9889 | 99.7% | 49,276 | 39:13:13 | 89.03 GiB |
| ProteinOrtho 6.3.6 | default ProteinOrtho mode; DIAMOND | 1.0000 | 0.8700 | 0.9305 | 96.0% | 37,201 | 0:09:29 | 2.57 GiB |
| SonicParanoid 2.0.9 | default mode; DIAMOND very-sensitive | 1.0000 | 0.7952 | 0.8859 | 87.8% | 19,853 | 1:14:39 | 7.90 GiB |
| OrthoHMM 0.5.0 | high sensitivity + inferred phylogeny (satellite_v2) | 1.0000 | 0.7733 | 0.8721 | 100.0% | 152,748 | 2:17:05 | 10.98 GiB |
| FastOMA 0.3.5 | final orthologous groups; supplied species tree | 1.0000 | 0.7569 | 0.8617 | 89.4% | 28,019 | 1:48:43 | 0.59 GiB |
| OrthoFinder 3.1.5 | full pipeline; root HOGs | 1.0000 | 0.7484 | 0.8561 | 88.1% | 67,131 | 1:45:32 | 3.99 GiB |
| OrthoHMM 0.5.0 | high sensitivity | 1.0000 | 0.7040 | 0.8263 | 100.0% | 164,826 | 1:49:16 | 11.12 GiB |

## Interpretation

OrthoHMM phylogeny improves F-score by +0.0458 over high sensitivity alone and is +0.0160 versus OrthoFinder's full root-HOG output. OrthoFinder's sequence-only checkpoint remains higher by 0.1174.

OrthoMCL is 0.0006 below OrthoFinder's sequence-only checkpoint, but its legacy BLAST stage makes it much slower on this dataset.

The benchmark is deliberately narrow: it scores gene-pair recovery only within conserved BUSCO families. Precision is often 1.0 because splitting a BUSCO family reduces recall while creating no cross-family false-positive pairs. This result should be interpreted alongside QfO and OrthoBench.

## Provenance Notes

- OrthoFinder sequence-only time is derived from the matching full 3.1.5 run at its MCL checkpoint; it is not a separate timed invocation.
- ProteinOrtho and SonicParanoid inference outputs were retained because they already used this exact biological input; they were rescored with the common evaluator for this report.
- Root-HOG outputs are reported for phylogenetic pipelines, while flat orthogroups are reported for sequence-only methods.
- OrthoMCL used its native legacy BLAST and inference rules. Its serial BioPerl BLAST-to-BPO conversion was replaced by a byte-compatible streaming converter validated against 6,417,790 native records, and independent species-pair calculations were process-parallelized before native matrix construction and MCL. Its wall time sums the separately measured BLAST, conversion, indexing, and accepted downstream stages. Failed validation attempts are excluded; the accepted 66-pair stage was rerun from zero checkpoints.
- OrthoMCL peak memory is the sampled sum of worker RSS values. Forked workers share copy-on-write pages, so 89.03 GiB is a conservative accounting value rather than physical unique memory.
- External-tool GNU time RSS values may omit memory held by container descendants and are not directly comparable to OrthoHMM's sampled process-tree RSS.
- Historical OrthoFinder 2.5.5 runs are excluded from the parity table; OrthoFinder 3.1.5 is the retained comparator.

Machine-readable details, checksums, job IDs, and source revisions are in `three_kingdoms_parity_20260907.json`.
