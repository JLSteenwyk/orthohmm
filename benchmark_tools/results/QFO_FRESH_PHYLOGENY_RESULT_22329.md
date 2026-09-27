# Fresh Retained-Order QfO Phylogeny

## Validation

Native job 22329 completed successfully in 01:53:08 on 32 allocated CPUs.
Corrected independent readback 22332 completed successfully in 00:22:44 on
2 allocated CPUs. Failed readback 22330 remains retained and documented in
[the inventory failure diagnosis](QFO_FRESH_READBACK_FAILURE_22330.md).

The [machine-readable comparison](qfo_fresh_phylogeny_result_22329.json) is
bound to the original native plan and all eight historical/fresh scientific
reports. The [receipt verification](qfo_fresh_phylogeny_readback_verified_22332.json)
rechecks the successful execution receipt, result, admission, reports, source,
and all 745 corrected-audit plan records (including 738 Python source files).
Both arms pass structure, sequence content, event/pair semantics and hierarchy
checks, with the full universe of 984,137 genes and 78 species.

## Historical Versus Fresh Retained Order

| Comparison | Result |
| --- | ---: |
| Identical root groups, label invariant | 366,068 |
| Historical-only / fresh-only root groups | 0 / 0 |
| Genes in changed root groups | 0 |
| Changed source-family partitions | 0 |
| Shared native ortholog pairs | 5,959,560 |
| Historical-only / fresh-only pairs | 0 / 0 |
| Changed species/confidence annotations | 0 |
| Species-tree bytes identical | Yes |

The fresh arm rebuilt trees without checkpoint reuse: 24,264 reconciled and
327,475 bypassed families, from 351,739 candidates, using 26 species-tree
marker families. Both arms report identical membership filtering: 40,169
constraints, 28,835 supported, 11,334 detached, 26,107 detached genes,
191,368 pairs removed and 11,186 root groups added.

## Scope And Next Step

This establishes reproducibility of the retained-order downstream predictions
under the recorded environment and inputs. It is not a new search run,
independent generalization test, accuracy improvement or dedicated-host timing
result. No historical scores have been replaced. In particular, equality for
retained order does not establish equality for the canonical-order candidate
arm, which changes two candidate groups and their constraints.

Next bind a tested reuse gate to successful audit 22332 while preserving the
audited source identities, then freeze and submit the canonical-order run.
Only input-identical raw-tree caches are eligible; species-tree inference,
rooting, reconciliation and membership filtering must run again. Compare all
three arms before determining whether any evaluation needs to be repeated.
The wider publication requirements remain open.
