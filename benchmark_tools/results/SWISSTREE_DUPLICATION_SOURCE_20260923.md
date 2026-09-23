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

## Complete Model-Source Acquisition

Subsequently acquired the uniquely linked current model tree from every one
of the 18 retained family pages. [Acquisition receipt](swisstree_model_acquisition_20260923.json)
retains URLs and hashes for all pages/trees; raw files remain under
`benchmarks/work/swisstree_model_sources_20260923`, outside Git. The downloader
does not execute page scripts, and refuses absent/ambiguous model links.
All 17 acquisition/inventory/inspection tests passed.

[Format and mapping feasibility](swisstree_model_feasibility_20260923.json)
retains all outcomes. POP, VATB, SUMF and APP fail the missing-or-duplicate
terminal-label check after Newick parsing; no deduplication or fabricated
accession mapping is applied. The other 14 parse successfully, yielding 366
exact-name/suffix mapping candidates across the 563-protein benchmark universe.
Those candidates are not admitted aliases. Only CASP and Clusterin phyloXML
files declare themselves rooted; NHX rooting semantics are not interpreted.

Explicit phyloXML duplication sums differ from the summary for CITE (3 versus
6), PSEN (8 versus 7), Clusterin (5 versus 4) and MAPT (1 versus 2).
Several trees contain no explicit duplication-count elements; their sum of
zero is not a biological claim of zero duplication. NHX event semantics are
uninterpreted and remain null. Missing implicit/alternative annotations,
source-version differences and sampling are possibilities, not established
explanations. Do not replace one source with the other to obtain preferred
counts. These discrepancies require semantic and benchmark-correspondence
review before freezing any ancestral-event feature or joining predictions.

## Multiline NHX Diagnostic

The four terminal-label failures are now localized: each source contains two
complete physical-line trees without semicolon terminators. The initial
Bio.Phylo whole-file parse combines these into one apparent tree with repeated
leaves. This is not evidence of biological duplicate identifiers. The original
failed feasibility report remains preserved rather than overwritten.

[Separate-candidate diagnostic](swisstree_multiline_model_diagnostic_20260923.json)
uses DendroPy 5.0.8 and appends a terminator only to each parser input, never
to the retained source bytes. Both candidates in each file have identical
leaf sets, stored clade sets and D=Y descendant sets. POP has 45 leaves and
two D=Y clades; VATB 55 and five; SUMF 35 and one; APP 45 and none. The APP
absence remains unknown event annotation, not a zero-duplication claim.
Support/color annotations need not be identical. Six focused diagnostic and
format tests pass. No candidate is chosen, and curated biological rooting,
accession mapping and benchmark pair-label correspondence remain unvalidated.
