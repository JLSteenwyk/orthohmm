# Original SwissTree Trees Retained in QfO

Follow-up: [native relation reconstruction](SWISS_RETAINED_MAPPING_20260923.md)
matches every family. It also identifies nine wrapped BAR annotations not
counted by this inventory's literal colon-token rule. The BAR value below
is therefore not a complete native duplication count; do not use these
inventory counts as outcome strata.

The frozen 2020 `ReconciledTrees_SwissTrees.drw` contains a `Tree`
argument in every `RecTreeCase`, not only the numeric relation table.
The pinned Darwin container reads these objects directly. The native
exporter traverses each binary tree in preorder; the Python reader checks
family coverage and traversal completeness without parsing Darwin tree syntax.

[Machine-readable inventory](swiss_retained_tree_inventory_v3_20260923.json)
records the reference/container/helper identities, native-output digest,
annotation counts and duplicate labels. Full trees and native stdout are
not redistributed. Original reference bytes remain unchanged.

| Family | Mapped proteins | Tree leaves | Explicit D=Y nodes |
|---|---:|---:|---:|
| APP | 27 | 125 | 3 |
| ASTER | 32 | 96 | 2 |
| BAMBI | 14 | 44 | 2 |
| BAR | 53 | 205 | 11 |
| CASP | 63 | 172 | 62 |
| CITE | 28 | 86 | 6 |
| Clusterin | 17 | 73 | 5 |
| GH14 | 15 | 159 | 39 |
| HOX | 56 | 90 | 25 |
| MAPT | 13 | 66 | 2 |
| NOX | 50 | 161 | 10 |
| POP | 34 | 105 | 3 |
| PSEN | 35 | 127 | 8 |
| RPS | 12 | 51 | 4 |
| SERC | 43 | 112 | 16 |
| SUMF | 12 | 24 | 1 |
| TRFE | 31 | 121 | 5 |
| VATB | 28 | 49 | 7 |

There are 563 mapped-protein memberships and 1,866 retained tree leaves.
POP repeats four labels twice each: ENSAMEG00000011433,
ENSECAG00000021367, ENSMGAG00000013529 and ENSTGUG00000012218.
The native generator explicitly handles intersecting child mappings;
therefore simply requiring unique labels or silently dropping duplicates
would not reproduce its semantics.

`AddReconciledTree.drw` defaults an internal node's event to S, then
examines its annotations. Missing annotations are consequently not explicit
speciation evidence. GH14 embeds D=Y in compound `:GN=...:D=Y` annotations;
token matching is required. An initial local inventory counted only exact
annotation strings and was superseded, not used scientifically.

These counts describe all retained tree taxa, not events restricted to the
benchmark proteins. They differ from current SwissTree downloads and do
not yet authorize duplication-history outcome strata. Next: recover the
native leaf mapping, preserve duplicate/intersection handling, reconstruct
relations against the retained relation table, and freeze feature definitions
before joining outcomes. No scores, thresholds or uncertainty estimates changed.

Reproduce from the repository root:

```bash
python -m benchmark_tools.inventory_swiss_retained_trees --repo . --output /tmp/swiss-retained-inventory.json
python -m pytest -q tests/test_inventory_swiss_retained_trees.py
```

The output path must not exist. Fifteen tests cover traversal failures,
duplicate preservation, missing annotations and compound NHX tokens.
