# Recovery Panel Coverage Gate

`benchmark_tools.verify_blast_recovery_panel` implements the complete replay
coverage check required by the
[recovery protocol](QFO_BLAST_RECOVERY_PROTOCOL_20260923.md).
It is prepared before all batches finish and has not admitted a live panel.

Before writing a report it requires all twenty native tasks in array 22103
and all twenty independent admission tasks in array 22105 to have completed
0:0 with their expected node and CPU allocations. It verifies the frozen
admission executor, manifest digests, per-batch report provenance and exact
ordered FASTA membership. Concatenated batch query identities must equal the
98,913-query replay FASTA without duplicates or omissions. Query block order,
coverage and diagnostics are checked again from the admitted records. All
distinct input/output records are hash-checked before and after aggregation;
conflicting records for a shared path are rejected.

The report retains failed queries separately from no-hit queries without a
logged failure. It cannot admit the interrupted prefix, merge BLAST tables,
release the old downstream chain or authorize scientific scoring. It relies
on the independent per-batch row audits rather than claiming a second native
search or wholly independent row-parser implementation.

```bash
python -B -m benchmark_tools.verify_blast_recovery_panel \
  --root "$PWD" \
  --output benchmarks/work/qfo_blast_recovery_panel_admission.json
```

Run only after all admission jobs complete; the output must not already
exist. No scheduler job for this aggregate gate has been submitted yet.
Forty-two focused panel/batch tests pass, including rejection of missing,
reordered, duplicate, unadmitted and prematurely promoted results, failed
queries with outgoing partial output, changed coverage and conflicting hashes.
Live whole-panel execution remains an explicit outstanding validation step.
