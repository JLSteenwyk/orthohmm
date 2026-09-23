# Mapped-Tree Duplication Features

The [protocol](SWISS_DUPLICATION_FEATURE_PROTOCOL_20260923.md) was committed
and pushed as **b288584 before feature extraction or outcome stratification**.
The [validated feature report](swiss_duplication_features_v2_20260923.json)
has SHA256 `97b0c4755d6a9df258d5c3f60fc0d5d25f1e5c09c42216c754a245a67d1942ec`.

Fresh native tree traversal and Python reconstruction again reproduce all
18 family mapping diagnostics and 10,765 reference pair labels. A separate
Darwin recursion, using the admitted leaf-to-entry assignments, agrees with
Python on mapped-member, duplication, explicit-speciation, default-speciation
and overlap counts for every family. This is a second traversal/counting
implementation, not independent biological validation or a second identifier
mapping source. No prediction results were loaded by the feature extractor.

There are **546 informative nodes: 128 explicit duplication annotations,
418 default-S nodes and zero explicit speciation annotations**. HOX contains
56 informative nodes for 56 unique mapped proteins because of its alias
collision. This is why the prespecified denominator is informative nodes,
not unique proteins minus one. Events with no mapped proteins on one side
do not count toward the feature.

The exact median duplication fraction is **7/48** (approximately 0.145833).
The frozen tie-preserving rule yields:

- Lower, at or below median: APP, ASTER, BAMBI, BAR, CITE, Clusterin, NOX,
  POP, SUMF (9 families).
- Upper, strictly above median: CASP, GH14, HOX, MAPT, PSEN, RPS, SERC,
  TRFE, VATB (9 families).
- Missing: none. Keep the empty missing bin in later score tables.

Fifty-five focused tests pass, including masked/unmapped clades, duplicate
aliases, informative-node denominators, wrapped annotations, exact rational
medians, ties, missingness and native-count disagreement. The first native
cross-check attempt failed closed because of Darwin table initialization/
indexing; the corrected full execution agrees, and no failed result was
admitted. A Python-only preliminary feature report remains local, superseded
by the dual-implementation report linked above.

Reproduce with a fresh output path:

```bash
python -m benchmark_tools.extract_swiss_duplication_features --repo . --output /tmp/swiss-duplication-features.json
python -m pytest -q tests/test_extract_swiss_duplication_features.py tests/test_audit_swiss_retained_mapping.py tests/test_inventory_swiss_retained_trees.py
```

The result is ready for the prespecified descriptive score join, not a
subgroup accuracy conclusion. The feature is reference-derived and related
to the reference's own pair labels. Default-S nodes are not explicit
speciation observations; missing events and alias collisions limit biological
interpretation. Development exposure and possible size/composition confounding
remain. No new scores, confidence intervals or superiority claims here.
