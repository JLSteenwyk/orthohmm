# Recovered Candidate Content Audit

`benchmark_tools/audit_blast_recovery_candidate.py` connects the existing
whole-table numerical audit to recovery query-disposition reconciliation.
It is implemented and fixture-tested, not executed on a production candidate.
All twenty replay batches and the guarded merge must finish first.

The auditor records input/helper identities, checks a contiguous original
prefix inventory with exactly one excluded terminal boundary query, requires
that boundary query in the replay universe, and derives expected per-query
row counts and SHA256 values. A fresh pass over the complete candidate checks
original FASTA query order, line completeness and exact block bytes. A separate
pass validates every HSP's identifiers, numbers, protein coordinates, alignment
accounting and contiguous query/subject blocks. The disposition check requires
the exact original query partition, no-hit identities, failed-query identities,
row totals and zero rows above the frozen cutoff. Input/helper identities are
rechecked before writing a fresh report.

This is deliberately not whole-search admission. It does not authenticate the
supplied batch reports' native execution, independently grant permission for
historical prefix reuse, or validate formatted-database sequence parity.
Those provenance and database checks belong to the still-required final
admission wrapper. All scientific/downstream authorization flags remain false.
No held dependency is released by this command.

Example interface, using verified paths and one `--batch` per admitted report:

```bash
python -m benchmark_tools.audit_blast_recovery_candidate \
  --blast /verified/merge/table/all.blast.candidate \
  --fasta /verified/original/all.fa \
  --log /verified/merge/selected.blast.log \
  --prefix /verified/original/query_blocks.jsonl \
  --batch /verified/batch_00_admission.json \
  --batch /verified/batch_01_admission.json \
  --output /fresh/candidate_content_audit.json
```

The illustrative two-report command is not the production twenty-batch
inventory. Production scheduling must also account for two complete scans,
before/after hashing, and in-memory query-block inventories; no runtime or
peak-memory claim has been established for this new audit.

Eighty focused recovery tests pass, including ten new integrated content tests
for correct output and corrupted bytes, order, completeness, numerical fields,
cutoff, failure records and boundary inventories. These fixture tests do not
replace the future full-scale audit.
