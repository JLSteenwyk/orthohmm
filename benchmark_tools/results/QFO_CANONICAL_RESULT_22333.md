# Canonical-Order QfO Phylogeny Result

## Validation

Native 22333 completed 0:0 in 00:12:16 with 32 allocated CPUs; independent
audit 22334 completed 0:0 in 00:14:40 with 2 allocated CPUs. The
[receipt verification](qfo_canonical_readback_verified_22334.json) rechecks
the final execution binding, all 740 audit-plan records, scientific report
identities/statuses and compared artifact identities. The canonical structure,
sequence, event/pair and hierarchy readers validate 984,137 genes in 78 species.
See the [complete three-way result](qfo_canonical_phylogeny_result_22333.json).

## Prespecified Contrasts

| Contrast (left versus right) | Left-only pairs | Right-only pairs | Shared-pair annotation changes | Identical root groups | Genes in changed root groups |
| --- | ---: | ---: | ---: | ---: | ---: |
| Historical versus fresh retained | 0 | 0 | 0 | 366,068 | 0 |
| Historical versus canonical | 110 | 85 | 7 | 366,067 | 44 |
| Fresh retained versus canonical | 110 | 85 | 7 | 366,067 | 44 |

Historical and fresh retained each have 5,959,560 pairs and 366,068 root groups.
Canonical has 5,959,535 pairs and 366,071 root groups. Each canonical contrast
shares 5,959,450 pairs; one old root group is replaced by four new groups.
The 44 genes are the union of genes in changed root groups, not a count of
genes individually reassigned. Source-family differences are restricted to
`Family0014602` and `Family0015042`. These IDs also reflect the previously
documented candidate-family reassignment, not two independent biological families.

All three species trees are byte-identical and have identical rooted topology
and taxa. Thus the observed differences are not explained by a changed species
tree. The experiment changes candidate membership and constraints together;
it does not isolate their individual causal contributions.

Canonical reused 24,263 input-identical raw-tree checkpoints, with zero remapped
checkpoints and no species-tree cache hit; it reconciled 24,264 families and
bypassed 327,475. Compared with retained order, membership filtering reports
three fewer supported constraints, three more detached constraints/genes,
24 more removed pairs and three more added root groups. These aggregate
differences are descriptive, not proof of improved biological accuracy.

## Interpretation And Next Actions

Retained-order reproducibility is established, but canonical ordering is not
prediction-equivalent on QfO. Preserve both arms and the negative equivalence
result. Do not transfer the old six endpoint scores to canonical predictions or
choose ordering based on which scores better. Next freeze conversion and scoring
against the same retained QfO references, tool versions and endpoint definitions,
then evaluate the validated canonical predictions. The unchanged species tree
does not justify skipping that evaluation.

GNU time reports 735.57 seconds wall, 1,160.76 user seconds, 118.01 system seconds
and 8,303,980 KiB maximum process RSS. This cached downstream run on a shared host
is not an end-to-end runtime, dedicated timing measurement, process-tree memory
peak or speed comparison with fresh inference. No historical accuracy score has
been replaced. Independent generalization, uncertainty, controlled timing and
the remaining publication requirements are not satisfied by this experiment.
