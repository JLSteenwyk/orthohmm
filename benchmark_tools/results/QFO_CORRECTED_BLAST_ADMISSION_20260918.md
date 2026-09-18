# Corrected OrthoMCL Search Admission

## Frozen Job

Validator revision `ae7de310cf7fbcba59107f2b054fc186f206eec4` is frozen in
`benchmarks/work/publication_qfo_corrected_blast_admission_v1`.
Job 21746 is queued `afterany:21713` on bizon with 2 CPUs, 64 GiB RAM,
24 hours and automatic requeue disabled. It rejects unsuccessful native
execution despite being scheduled after any terminal outcome.

The native search remains frozen at
`b0c7128d464e2abe5bf1fd3646c0378625bbb0de`; its legacy engine, masking,
commands and 180-CPU/900-GiB allocation are unchanged. The validator is not
a search restart and does not touch unfinished search outputs.

## Checks

`admit_qfo_corrected_blast.py` requires uniquely completed scheduler evidence
with exit `0:0`, the expected node/allocation, a clean pinned search executor,
the frozen input/runtime manifests and exact execution provenance. It checks
native stage exit codes, finite ordered timestamps, log/timing paths, final
database inventory and final BLAST table identity. A remaining partial table
alongside the final output is rejected.

It then extracts the complete formatted database using the frozen legacy
runtime and requires exact input/header/order parity for 984,137 sequences.
No case conversion, terminal-stop removal or other residue normalization is
silently accepted. Differences retain a component report and fail the combined
gate for review.

The complete BLAST table and diagnostic log are then audited for structural
validity and query/subject/self-hit coverage. Logged failed queries remain
explicit in the admitted evidence; `search_admitted: true` does not assert
that every query succeeded. Missing hits without logged errors are not
automatically classified as failures. Original-release failure counts are
not reused.

All bound artifacts and generated audit outputs are rechecked before issuing
`corrected_orthomcl_search_evidence_verified`. The future report directory is
`benchmarks/work/qfo_corrected_blast_admission_20260918`, containing:

- `report.json`: combined provenance and search-evidence decision.
- `database/report.json`, `database/database.fasta`, `database/fastacmd.log`:
  extraction and sequence-parity evidence.
- `table.json`: complete hit-table structure and diagnostic query coverage.

Errors after component validation begins produce a retained failed combined
report. Earlier scheduler/provenance failures appear in the validator job log
without fabricating a completed audit. Existing output directories are never
implicitly resumed or overwritten.

## Verification And Limits

The complete unit suite passed with installed legacy-engine smoke tests
enabled: 3,628 passed. The 33 new driver tests cover scheduler/provenance
drift, native stages, failure preservation, database differences, input drift,
and no access to unfinished jobs. Driver orchestration tests mock component
audits; real component smoke evidence is separately retained in
`ORTHOMCL_DATABASE_AUDIT_20260918.md` and
`ORTHOMCL_SEARCH_TABLE_AUDIT_20260918.md`. Bash syntax validation passed.

At submission, native search 21713 was resource-pending and validation 21746
was dependency-pending. No complete corrected search has yet been admitted.
BPO/index validation, native final groups, pair conversion, reference impact
and accuracy scoring remain separate unfinished stages. This report does not
authorize downstream execution or claim matched timing or publication readiness.
