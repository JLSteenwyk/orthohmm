# SwissTrees Frozen Mapping Coverage Audit

Following the exact-accession inventory, audited all563reference accessions
against the frozen2020QfO numeric protein mapping and every accession in
the78frozen input FASTAs. No similarity search, current external mapping,
fuzzy name match or best-scoring candidate is used.

```sh
python benchmark_tools/resolve_swiss_sequence_aliases.py \
  --inventory benchmark_tools/results/swiss_sequence_inventory_20260917.json \
  --prepared benchmark_tools/results/qfo_factorial_prepared_20260917.json \
  --mapping qfo_benchmark/benchmark-webservice/reference_data/2020/mapping.json.gz \
  --output benchmark_tools/results/swiss_sequence_alias_audit_20260917.json
```

Report SHA-256:
`b69f74e35aa3fa3af99a33b3b27f830d79df2abd3935b84e911f981186ac33ec`.
The mapping itself is pinned to
`1c10f6ce5e53ebc3148dde02d16268b225c3c8817c952d1274daa41acbf9eb4d`.

## Outcome

- 549 reference proteins match unique exact input accessions.
- Zero resolve to a different input accession through the same numeric ID.
- Zero have ambiguous input matches.
- Fourteen numeric reference IDs have no matching input accession.

The14missing identities and their numeric IDs are retained in the report;
their accession list is also in `SWISS_SEQUENCE_INVENTORY_20260917.md`.
All14map to the retained domain annotation file `UP000008143_8364.json`,
as checked against `swiss_domain_annotation_inventory_20260917.json`.
The corresponding input FASTA identifies Xenopus tropicalis, taxon8364.
This localizes the missing identity coverage; it does not establish why
the input/reference resources differ or exclude an unmapped homolog.

All549previous exact-match descriptors are unchanged. The missing14cannot
be filled by the frozen numeric alias mapping. Their sequence descriptors
remain missing; reference families and benchmark score denominators are
not silently reduced. No inference, reference truth or scoring settings change.

## Consequences And Outstanding Checks

This is an input-coverage limitation for the frozen factorial FASTAs.
It is not yet a quantified correction to retained accuracy, nor proof that
all historical comparator inputs share the same limitation. Audit affected
reference pairs and historical input identities before attributing score
differences or asserting a shared recall ceiling. Preserve the original
full-reference endpoints; any available-input sensitivity analysis must be
separately declared and must not replace the frozen primary result.

No positive explicit fragment labels were recovered. The absence of labels
still does not establish sequence completeness. Subsequent composition
stratification must prespecify handling of partially covered families.

The resolver rejects duplicate input accessions and nonunique reference
numeric identities, reports ambiguous matches without choosing one, and
checks source hashes before/after reading. Tests include exact/alias/missing/
ambiguous cases, invalid IDs, a file-level alias fixture, preserved exact
descriptors, and changed reference membership.
