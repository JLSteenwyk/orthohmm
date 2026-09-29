# Corrected VGNC Cross-Block Components

The [eight-method table](CORRECTED_VGNC_COMPONENT_TABLE_20260928.md) and
[machine-readable result](corrected_vgnc_components_20260928.json) describe
connected components induced by scored cross-block false positives. This
extends the historical four-stage dependency description to the frozen
corrected comparison: all eight individual graphs, all 28 pairwise unions,
and the all-method union are retained. No confidence interval is calculated.

## Construction

Start with the 16,844 reference-defined overlap blocks and their 23,934
asserted pairs. Read the previously validated sparse TP/FP/FN count tables.
Each nonzero off-diagonal cell adds one undirected link between two blocks,
regardless of its FP multiplicity. Within-block false positives stay in the
count totals but add no graph link. Every method must preserve each block's
TP+FN count and the previously audited overall count triple. Duplicate cells,
unknown blocks, negative counts or cross-block truth are rejected.

Reference blocks without prediction links remain singleton components. For
paired comparisons, take the union of the two methods' links before finding
components; do not select the more favorable method's graph. The all-method
union uses every retained comparator, including sequence-only OrthoFinder.

## Observations

- Phylogenetic OrthoHMM has eight cross-block links, connecting 16 blocks in
  eight two-block components. Its 13 FP observations are not eight FP pairs:
  multiple scored pairs can occupy the same block-pair cell.
- Full OrthoFinder has 38 links; its largest component contains four blocks.
- Their paired union has 41 distinct links, 65 linked blocks and 16,808 total
  components. Its largest component has four blocks and 24 asserted pairs;
  176 asserted pairs lie in nonsingleton components overall.
- The all-method union has 136,107 distinct links and 6,615 components,
  including 4,187 isolated blocks. Its largest component has 370 blocks and
  470 asserted pairs. Nonsingleton components contain 12,657 blocks and
  17,987 asserted pairs. Thus pooling comparators changes the proposed units
  substantially, even though the reference itself has not changed.

These are dependencies visible in the scored pair tables, not a complete
biological dependence graph. The component with the most blocks need not
contain the most truth pairs; the result retains both size and truth mass.

## Statistical Boundary

Connected-component resampling is not automatically justified by this graph.
The components depend on observed method predictions, and changing the set
of comparators changes them. Shared evolutionary history or correlated errors
can connect otherwise disjoint components without any scored cross-link.
Conversely, an FP link alone does not specify a stochastic dependence model.
The rare-error and shared-clade coverage failures remain unresolved. This
audit supplies no independent sampling law, bootstrap validity or native CI.

Twenty-seven focused tests pass across this inspector and the preceding block
mapping/dependency tools. A separate standard-library breadth-first traversal
reproduced the size histograms, link/component counts, largest-component sizes
and reference-pair masses for [all 37 graphs](corrected_vgnc_components_crosscheck_20260928.json).
All 40 source/input record entries were rechecked around that cross-check. No native prediction,
score, reference category, default or uncertainty claim changed.

```bash
python -m benchmark_tools.audit_corrected_vgnc_components \
  --source benchmark_tools/results/corrected_vgnc_blocks_20260926.json \
  --output /fresh/path/corrected-vgnc-components.json
```
