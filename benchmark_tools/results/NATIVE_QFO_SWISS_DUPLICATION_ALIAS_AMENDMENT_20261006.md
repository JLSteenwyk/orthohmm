# Native Duplication-Strata Identifier Join Correction

Initial real exporter exits1 before output/scoring with "Native and retained
mapped-tree memberships differ". The incorrect code compared retained
`mapped_labels` keys (original tree leaves, including Ensembl aliases) directly
with native raw UniProt accessions. This is an implementation error, not
evidence of a scientific reference mismatch; investigate entry IDs before
admitting the join. Preserve the original time receipt and source1b893235.

Correct join: directly check/read the original retained `mapping.json.gz`
(SHA256 `1c10f6ce5e53ebc3148dde02d16268b225c3c8817c952d1274daa41acbf9eb4d`)
and map each native accession to its positive integer reference-entry ID.
Require distinct IDs for distinct native genes, exact equality with the set
of retained `mapped_labels` values, exact prior mapping status and mapped-member
counts for every family. Retained leaf aliases can share an entry ID (HOX).
Do not strip labels, guess aliases or use count-only equivalence. Both readers
recheck the identifier resource; all other original traversal/acquisition
evidence remains inherited. No tree traversal or prediction re-inference.

Add regression tests for Ensembl aliases, duplicate retained leaf aliases,
unknown/wrong/invalid native entry IDs, native collision and missing family.
Commit correction/tests before recomputing native bin outcomes. Preserve
September23 fraction definition/median/ties/bin membership and the native
statistic, endpoints and inference settings. Output path may still be used
because the failed run creates no output; no timing inference retry occurs.

Initial test invocation in scientificPython3.10 also fails before collecting
tests because pytest is absent. Do not install there. Forty tests pass0.59s
and121 joined tests pass1.88s in retained testPython3.12 before the first real
run; keep those receipts as pre-correction history, not final tested source.
