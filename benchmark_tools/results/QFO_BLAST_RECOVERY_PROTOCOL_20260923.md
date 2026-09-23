# Corrected OrthoMCL Checkpointed Search Recovery

## Evidence And Scope

Diagnostic 22055 completed 0:0. Its hash-bound comparison reports compatible
combined/isolated HSP multisets and diagnostics for all five queries. Three
earlier complete query blocks match exactly (20, 20 and 126 HSPs). The known
failed query again has no HSPs and identical failures. The interrupted final
query has 73 replay HSPs versus 47 retained complete rows; those 47 are a
multiset subset, not a complete result. The entire old final block is excluded.

Comparison SHA-256:
9d49960fd22aaa1b538dcae728b06d0a67802d8d3d1f4296b31c459a98f702e1.
Partition manifest SHA-256:
c37381bef1cdb406849f62accde679480a4abe76193cc0fdc65a90e16c8d6a34.
The old 33,141,805,056-byte partial output remains immutable. Neither five-query
agreement nor this plan alone admits the retained prefix or final search.

## Replay Execution

Replay all 98,913 unresolved input queries in original order, including all
98,912 absent queries and the interrupted query. Split the exact-byte replay
FASTA into 20 consecutive batches: 19 of 5,000 queries, then 3,913. Do not
select batches, reorder by score or discard known failed/no-hit queries.
Every batch searches the unchanged complete 984,137-protein formatted
database. Never format a smaller database from these query files.

Preserve the installed legacy BLAST 2.2.13 binary, runtime snapshot, original
working directory, clean environment, and all search flags. Only -i and -o
change. In particular -p blastp, -e 1e-5, -m 8, -a 180, -v 1000 and -b 1000
remain frozen. Use a serial scheduler array, 180 CPUs and 900 GiB per task,
24-hour limit, no automatic requeue. These are recovery costs on a shared
workstation, not controlled comparative timings. No DGX work is authorized.

Freeze and test the executor before submission. Each batch uses a fresh
directory and records command, input/database/runtime/source identities,
scheduler identity, start/end times, exit code, diagnostics and /usr/bin/time
output. Refuse pre-existing execution artifacts; no implicit restart. Flush
and fsync complete native outputs and logs before an atomic completed-status
record, and fsync its directory. A completed native process is still pending
query/output admission. Sequence-specific failures remain explicit even with
exit zero. Do not silently retry or replace a failed batch; preserve it and
review a separately recorded recovery if necessary.

## Consolidation Gate

Before any merge, independently require all 20 scheduler/native completions,
exact query/database/source identities and disjoint full replay coverage.
Validate every batch's HSP structure, complete query blocks and diagnostics;
every emitted query must belong to that batch and every subject to the full
database. No-hit queries need completed execution evidence, not invented rows.
Preserve query failures and distinguish outgoing search failures from incoming
subject hits. Compare repeated diagnostics for the previously failed queries.

Retained-prefix reuse requires an explicit evidence review combining the
archived sequential-query source path, exhaustive row/block audit, unchanged
full-file hash, completed diagnostic agreement and query partition. Recheck
all prefix byte ranges and hashes. If this review does not support reuse,
search the remaining queries separately; do not lower validation standards.
Keep any residual source-to-binary/durability inference limitation explicit.

If reuse is supported, write a new table in original query order: complete
retained blocks only before byte 33,141,800,005, plus full replay blocks for
the replay set. No duplicated query blocks, no concatenation that ignores
original ordering, no old final-query rows, no merging failed partial batches.
Validate the entire resulting table and exhaustive query dispositions before
BPO conversion, indexing and clustering. Record combined provenance and
incremental costs, never relabel interrupted elapsed time as successful timing.

Held job 21746 and its downstream chain remain untouched until a separate
reviewed admission replaces their original unresumable search expectation.
Final OrthoMCL scores, failure-impact analysis and publication readiness are
not established by replay completion alone.
