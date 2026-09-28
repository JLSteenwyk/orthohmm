# TreeFam-7 A+B Archive Inspection

Downloaded the [TF7AB SQL archive](https://selectome.org/ftp/MySQL/selectome_04-TF7AB__mysql5.0.sql.zip)
listed by the [Selectome public archive](https://selectome.org/ftp/MySQL/).
Its 59,845,652 bytes match the published SHA256
`d936fa7fd1bbaa96904d787c705f54f2aeab4ad936f3ab4cd9ace41f0595d6e0`.
The sole SQL member is 371,301,900 bytes, SHA256
`c85aa1dd0d820deffc6d55d1abbc8972ba2dac78f495d4a132aa09e17862fbe9`.
ZIP integrity passes. The embedded news identifies a TreeFam7 A+B-based
Selectome release dated August 18, 2010. Raw files remain outside Git in
`benchmarks/work/treefam_selectome_search_20260928/`.

## New Material

The [streamed inventory](selectome_tf7ab_inventory_20260928.json) establishes:

| Item | Count |
| --- | ---: |
| Gene-table records | 1,287,270 |
| Gene taxonomy IDs | 67 |
| Parsed NHX subtrees | 9,850 |
| Distinct tree family accessions | 8,355 |
| Tree leaf occurrences | 168,719 |
| Earlier A-only subtree keys retained | 1,211 / 1,211 |
| Identical NHX strings among those keys | 85 / 1,211 |

All subtrees are labeled Euteleostomi. All gene-table primary indices and
sequence IDs are unique within the inspected table. The gene schema contains
sequence/transcript/gene identifiers, version numbers, display labels, taxonomy
and descriptions, not a sequence column or the original QfO mapping. No exact
duplicate leaf labels occur within any subtree. Counts of leaf occurrences are
not counts of distinct genes across all trees.

The 1,126 non-identical shared NHX strings have not been decomposed into
annotation versus topology changes. Matching family/taxon/subtree-number keys
and a common release label do not establish interchangeable tree contents.
No `treefam2reference` filename string appears in the complete SQL payload.

## Validation

SQLGlot 30.20.0 parses selected MySQL INSERT statements as literal data only.
The archive is streamed; no database is created and no downloaded SQL statement
is executed. Every gene row's width and identifiers are checked, and each tree
is parsed with preserved leaf spelling/case and required species annotations.
The original A-only inspector is unchanged.

All 13 archive/parser tests pass, including published checksum mismatch, wrong
earlier source, unexpected members, escaped strings and nonliteral statements.
An [independent Biopython NHX pass](selectome_tf7ab_tree_crosscheck_20260928.json)
agrees with DendroPy on tree/family/leaf totals and subtree-taxon counts, and
reproduces the whole SQL member digest. It shares the SQL literal parser and
does not independently validate the gene rows. The full main pass completed
successfully; no incomplete result was promoted.

With SQLGlot 30.20.0 and DendroPy available:

```bash
python -B -m benchmark_tools.inspect_selectome_tf7ab --directory benchmarks/work/treefam_selectome_search_20260928 --earlier benchmark_tools/results/selectome_tf7a_inventory_20260928.json --output /tmp/selectome_tf7ab_inventory.json
```

## Remaining Boundary

This recovers additional historical identifiers and derived vertebrate trees,
not an authenticated complete TreeFam-A collection or `treefam2reference.txt`.
No QfO pooled-reference reconstruction, family assignment, uncertainty estimate,
new score or reference substitution follows. Identifier and tree suitability
would require explicit validation; nothing here proves the old and new NHX
strings encode the same topology or event annotations. No person was contacted,
scientific runtime changed, unrelated job modified or DGX accessed.
