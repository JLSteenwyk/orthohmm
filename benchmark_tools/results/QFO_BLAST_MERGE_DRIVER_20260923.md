# Guarded Recovery Merge Driver

`run_blast_recovery_merge.py` connects the tested streaming merge component
to the [conditional prefix review](QFO_BLAST_PREFIX_REUSE_REVIEW_20260923.md).
It is prepared but **not submitted or executed on production BLAST bytes**.
The held whole-search admission and downstream OrthoMCL chain remain unchanged.

Before any copying, the driver requires completed 0:0 accounting with correct
node/CPU allocations for prefix recheck 22148, all twenty native batches and
all twenty batch validators. It pins the review, prefix recheck, query
partition, diagnostic comparison and merge/panel implementations. It checks
the frozen recheck executor, archived source/release identities and installed
binary equality with the archived executable without running downloaded code.

It rechecks recorded inputs and runs the full independent replay-panel gate
afresh. It verifies exhaustive query ordinals against the original FASTA,
the exact replay order, and consistent diagnostic categories/messages for
queries appearing in both the interrupted and replay logs. Original failed
queries cannot be treated as retained successful searches. Duplicate-message
multiplicity matters; log line numbers need not agree.

On success it writes a durable preflight record, interleaves hash-checked
blocks in original order and checks the final HSP-row total. The selected
diagnostic log contains unchanged original lines only for retained queries
and unchanged replay lines for replayed queries, in original query order.
It is not an unfiltered concatenation that double-counts repeated diagnostics.
The original logs remain preserved and hash-bound. Source identities are
checked again before durable terminal status.

The output is `table/all.blast.candidate` with status
`merged_candidate_pending_full_table_admission`. `search_admitted`,
`reuse_authorized` and `publication_ready` remain false. Exceptions retain
partial artifacts and an explicit failed status; no implicit resume or
overwrite is supported. An independent full numerical/alignment/diagnostic
audit and exhaustive query-disposition reconciliation are still required.

Forty-five scan/merge/driver tests pass, including an exact synthetic
candidate transaction, missing/failed/wrong-resource prerequisites, changed
review, diagnostic discrepancies, no-overwrite, corruption and fsync failure.
A real live prerequisite check rejected incomplete batch 22103_3 before
creating any merge directory or reading the production BLAST table. No
full-scale performance, complete preflight or real merge success is claimed.

Next: freeze a tested scheduler executor after the remaining integration
review, submit only with the complete-panel gates intact, then independently
validate its candidate before reconnecting conversion. Full native recovery
costs remain shared-host observations, not matched comparative timings.
