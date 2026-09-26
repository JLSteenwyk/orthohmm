# Hierarchy and Reconciliation Selection Readback

`audit_phylogeny_hierarchy.py` consumes the prior event/pair semantic report.
It independently checks that reconciliation was selected exactly for families
with at least three genes, at least two species and at least one multicopy
species. All other families must be bypassed. It reconstructs every hierarchy
row from validated internal node records, adding one unambiguous root row per
bypassed family. Table schema, order, identifiers, parents, mapped species-tree
nodes, events and full gene membership must agree exactly, with no extra or
missing rows.

## Results

Both the installed fixture and historical full OrthoBench p1_c1_r1 pass.
The historical result contains:

| Quantity | Verified count |
|---|---:|
| Reconciled families | 8,681 |
| Bypassed families | 45,764 |
| Hierarchy rows from internal nodes | 163,395 |
| Bypass root rows | 45,764 |
| Total hierarchy rows | 209,159 |

[Historical receipt](installed_ob_historical_hierarchy_20260926.json), SHA-256
`7a9ce8a72559bb3aebbe73197575be74bd3da5b479bc3967d2b927b59f565899`,
and [fixture receipt](installed_phylogeny_hierarchy_fixture_20260926.json).
All 100 related reader/runner tests pass, including selection-boundary cases,
incorrect reconciliation selection and malformed or incomplete hierarchy rows.

## Output Semantics

The hierarchy table represents the reconciled node hierarchy **before satellite
membership filtering**. It is not the final root-HOG partition and need not
encode its detached satellites as separate hierarchy nodes. The benchmarked
root groups and native pair table reflect membership filtering. This is the
frozen writer's verified behavior, not an inferred new meaning or a change to
scoring. Do not substitute one output for another.

These checks establish reconciliation/bypass selection and serialization
consistency, not biological correctness or correctness of upstream candidate
selection. They are conditional on the pinned inputs and the preceding event
readback. Job 22179 is still running and has not been admitted; full-dataset
installation reproducibility remains unproven until its results are checked.
The [reproduction guide](../PUBLICATION_REPRODUCTION.md) now lists the complete
post-completion readback sequence. No scores or scientific defaults changed.
