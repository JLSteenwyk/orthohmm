# Installed OrthoBench Run: Valid Outputs, Historical Disagreement

## Outcome

Job 22179 completed once, exit 0, with scheduler elapsed 2:43:10 and 32
allocated CPUs. All five prespecified readers passed: frozen provenance and
partition scoring; tree/pair structure; candidate/alignment/supermatrix sequence
content; independent event and membership-rule reconstruction; hierarchy.
The [compact receipt](installed_orthobench_readback_20260926.json) retains
reader identities, complete family score records, and hashes of full local
receipts. Large per-file inventories remain in the private run directory.
Its SHA-256 is
`79fc6ce6a169d486235c811c7d6eec9019127958c07e56d8ff982675cd4149e0`.

This is **not exact reproduction**. No historical scores, default parameters,
or comparator tables were replaced. No inference retry was performed.

| OrthoBench statistic | Historical P1/C1/R1 | Fresh installed run | Difference |
|---|---:|---:|---:|
| Weighted RefOG F1 (%) | 74.106074 | 73.821569 | -0.284504 pp |
| Precision (%) | 81.770454 | 81.103625 | -0.666829 pp |
| Recall (%) | 67.755336 | 67.739444 | -0.015893 pp |
| Root groups | 59,770 | 59,812 | +42 |

Both partitions cover all 251,378 input genes exactly once. They share 56,652
groups irrespective of labels; 3,118 historical and 3,160 fresh groups differ,
involving 26,433 genes. Sixty-two of 70 reference-family records are identical.
Changed records are RefOG007, 011, 021, 023, 027, 035, 053 and 058. Their exact
TP/FP/FN and split changes are retained in the receipt, including improvements
and deteriorations. These record differences are not additional independent
observations or a new selection endpoint.

## Validation Scope

The fresh run has 54,495 candidate families, 8,624 reconciled families and
45,871 bypassed families, with no checkpoint reuse. The structure reader parsed
25,873 trees. Sequence checks verified 8,624 gene alignments, 200 species-tree
alignments and the 179,564-column supermatrix. Independent reconciliation-rule
readback reproduced 965,660 native ortholog pairs (843,902 high and 121,758
medium confidence) and the final root partition. Hierarchy readback matched
all 209,176 rows, explicitly before satellite-membership filtering.

These are consistency checks conditional on saved trees and constraints. They
do not prove biological correctness, optimal rooting, or correctness of the
initial search. The audit environment shares the DendroPy parser family with
the implementation. All 85 focused tests for the five readers pass.

## Resources and Next Investigation

`time -v` reports 9,787 displayed wall seconds, 281,292.64 user seconds,
3,740.24 system seconds and 6,154,024 KiB maximum process RSS. These are
shared-host observations, not dedicated comparative timing or aggregate
concurrent-process peak memory. The earlier invalid `sstat` CPU observation
is excluded. No DGX use occurred.

The frozen scientific source was installed with changed numerical/dependency
versions and rebuilt MAFFT, as prespecified in the
[protocol](INSTALLED_ORTHOBENCH_PROTOCOL_20260926.md). This does not establish
the cause of disagreement. Next compare saved search, graph, initial/final
partitions, candidate families and species-tree inputs to locate the earliest
divergence before considering a targeted diagnostic. Do not attribute the
change solely to phylogenetic reconciliation or tune to recover the old score.
Full-dataset numerical equivalence remains unestablished; broader publication
requirements remain open.
