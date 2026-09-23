# SwissTree Duplication Source Feasibility

Retrieved the official [SwissTree family table](https://swisstree.sib.swiss/cgi-bin/swisst?page=gold_standard)
on 23 September 2026, without joining method predictions. A structured HTML
parser retained all 19 rows, names, published duplication counts and update
markers. The [machine-readable inventory](swisstree_duplication_source_20260923.json)
records the source SHA-256 and an explicit name crosswalk to the 18 retained
QfO families; ST012 (ant transformer) is excluded because it is absent from
that panel. Seven parser tests pass. Counts are source summaries, not newly
validated benchmark-specific ancestral histories.

Raw HTML: `benchmarks/work/swisstree_duplication_source_20260923/gold_standard.html`,
11,633 bytes, SHA-256
`928a6f09de5c6f9fcc32b02aea020ce722a8ada794bcf6361d63cf75eb7aadd0`.
The current site marks POP, VATB, SUMF and APP as updated. Its sampling and
release need not equal QfO 2020. Raw site files remain outside Git while
redistribution terms are unresolved; the page asserts SIB copyright.

## Tree Pilot

The official ST013 page links its
[Transferrin reference tree](https://swisstree.sib.swiss/ST/ST013/ST013_treemodel.phyloxml).
Downloaded it to
`benchmarks/work/swisstree_duplication_source_20260923/ST013_treemodel.phyloxml`;
SHA-256 `425a39c745cbf2a4a1967db961fa01ddfe20a366c3f00eff9e29977d9f658806`.
Bio.Phylo parsing finds 121 terminal leaves and five explicitly annotated
duplication events, agreeing with the site's ST013 summary. Splitting leaf
names at their final underscore matches all 31 retained TRFE accessions.
This is a feasibility check, not an admitted general alias-mapping rule.

The file declares `rooted="false"` despite displaying oriented clades and
event labels. Do not silently root the tree or interpret unannotated nodes
as speciations. Exact event provenance, rooting semantics, leaf uniqueness,
benchmark mapping and pair-label concordance remain to be established.
No duplication-stratified scores or bins have been generated.

## Next Checks

Acquire linked model trees for all retained families, keeping updated and
historical versions distinct. Validate leaf mappings against retained QfO
identities, then compare explicit event-derived relationships with the native
QfO reference table. Resolve rooting and unknown-event handling before any
ancestral-history feature. Freeze any usable feature definition before
joining outcomes; otherwise retain unavailable families as missing.

These curated annotations are independent of OrthoHMM/OrthoFinder predictions
but share the SwissTrees reference underlying this benchmark. They are not
independent validation data or experimental ground truth. Published event
totals also depend on tree size and taxon sampling. A source inventory does
not close the independent duplication-history analysis requirement.
