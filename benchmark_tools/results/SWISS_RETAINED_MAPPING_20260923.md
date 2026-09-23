# Retained SwissTree Mapping Reconstructed

All **18 families, 563 mapped memberships and 10,765 reference pairs**
reconstruct exactly from the retained tree objects and the frozen QfO
`mapping.json.gz`. Every family has zero missing/extra members, zero
missing/extra pairs and zero differing event labels. The
[machine-readable report](swiss_retained_mapping_v2_20260923.json) preserves
per-family comparisons, mapped labels, unmapped-leaf counts, intersection
observations, input/helper hashes and a digest of native stdout.

This uses Darwin itself to read the original tree and relation objects.
A separate Python traversal maps exact leaf identifiers first, then the
first matching underscore-delimited alias, as the retained generator does.
It drops overlapping identifiers from the left child, writes cross-child
relations, and permits later writes to replace earlier relations, matching
the generator's order. The reference relation table supplies the comparison,
not any tool's predictions. No new raw trees are redistributed.

The retained original `IDIndex.db` is unavailable at the reference path;
the test instead uses the supplied identifier mapping with SHA256
`1c10f6ce5e53ebc3148dde02d16268b225c3c8817c952d1274daa41acbf9eb4d`.
Exact agreement supports this mapping for these trees, not equivalence of
the entire missing index. The reconstructed relations are not independent
biological confirmation of the curated reference.

## Annotation and Duplicate Details

The generator recognizes event substrings inside annotations, including
BAR's `3[&&NHX:D=Y]` form. The earlier inventory's colon-token count excludes
nine wrapped BAR annotations: its value 11 is not the complete native
duplication count. That inventory remains a literal-token diagnostic and
must not be used as a duplication stratum. The reconstruction uses the
generator's event interpretation and reproduces all BAR relations exactly.

HOX contains one intersecting child mapping, entry 902573, even though its
literal leaf labels are unique. POP's repeated literal labels do not produce
a mapped-child intersection in this retained benchmark. Thus literal label
uniqueness alone is insufficient to validate the mapping. These differences
were retained rather than normalized away.

Missing annotations still default to S under the benchmark generator;
they are not explicit speciation observations. Event reconstruction does
not authorize a biological claim about all unannotated nodes. No accuracy
stratification or outcome-dependent threshold was selected here.

## Verification and Next Step

Thirty-nine focused tests cover native traversal, duplicate preservation,
mapping priority, wrapped annotations, event precedence, left/right overlap,
overwrite order, malformed reference exports and deliberately mismatched
members/pairs/events. Unsupported annotations fail closed. Reproduce:

```bash
python -m benchmark_tools.audit_swiss_retained_mapping --repo . --output /tmp/swiss-retained-mapping.json
python -m pytest -q tests/test_audit_swiss_retained_mapping.py tests/test_inventory_swiss_retained_trees.py
```

The output path must be new. The next analysis must freeze a duplication
feature on the benchmark-mapped tree, distinguish retained explicit events
from default-S nodes, and state duplicate-handling and tree-size normalization
before joining method outcomes. Existing development exposure remains;
this does not convert QfO into an independent test set.
