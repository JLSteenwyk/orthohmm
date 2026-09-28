# Shared Selectome Trees Preserve Compared Content

The preceding TF7AB inventory found only 85 byte-identical NHX strings among
1,211 family/taxon/subtree-number keys shared with the TF7A archive. A complete
structural comparison now resolves that discrepancy for the measured fields.

| Compared field | Shared records differing |
| --- | ---: |
| Case-sensitive leaf identities | 0 / 1,211 |
| Rooted topology (descendant clades, multiplicity retained) | 0 / 1,211 |
| Branch lengths (exact parsed values) | 0 / 1,211 |
| Duplication-event `D` annotations by clade | 0 / 1,211 |
| Species `S` annotations by clade | 0 / 1,211 |
| All parsed annotations by clade | 0 / 1,211 |
| Internal node labels by clade | 1,126 / 1,211 |

The 85 exact NHX matches have identical signatures in every measured field.
For the remaining 1,126, detected differences are confined to internal labels;
there is no detected topology, branch-length or annotation difference. Child
order is ignored, while case-sensitive leaf identity and repeated descendant
clades are retained. This comparison does not assert identity of every possible
Newick serialization detail.

The [complete comparison](selectome_shared_tree_comparison_20260928.json)
retains all keys and differing-field lists. Both archive hashes are pinned and
rechecked; their sole SQL members are read to EOF with ZIP CRC verification.
Only subtree INSERTs are parsed as data. No database statements are executed,
no tree is rerooted and no gene identifiers are mapped or case-folded.

## Verification

All 18 focused tests pass. Synthetic controls distinguish child-order invariance,
internal-label/support changes, branch-length changes, topology changes,
duplication-event changes and leaf case changes. Exact duplicate leaf labels
are rejected. Existing SQL and archive guard tests also pass. A second complete
execution reproduces the retained JSON after serialization.

```bash
python -B -m benchmark_tools.compare_selectome_trees --directory benchmarks/work/treefam_selectome_search_20260928 --output /tmp/selectome_shared_tree_comparison.json
```

Use the previously recorded SQLGlot 30.20.0 and DendroPy environment. This
analysis reads the two downloaded archives but does not rerun the expensive
full gene-table parse.

## Scope

This establishes consistency of compared content for shared derived vertebrate
subtrees. It does not authenticate complete original TreeFam-A 7 trees, recover
`treefam2reference.txt`, establish a mapping to QfO 2020 proteins, reconstruct
the pooled reference or prove biological tree correctness. No benchmark pairs
were assigned to families and no scores or uncertainty estimates changed.
Missing mapping and full-tree coverage remain the next source-suitability gaps.
