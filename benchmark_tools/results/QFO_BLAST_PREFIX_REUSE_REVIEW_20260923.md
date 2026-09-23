# Conditional Retained-Prefix Reuse Review

## Decision and Scope

The combined evidence supports reuse of the **885,224 complete query blocks
strictly before byte 33,141,800,005**, conditional on the executable checks
below. This is a documented recovery inference, not proof of historical
filesystem durability or a reproducible source-to-binary build. No complete
search, merged table, BPO input, orthogroups or final score is admitted here.

The proposed prefix contains 366,012,616 HSP rows. The entire final query
`tr|F1MRE1|F1MRE1_BOVIN` is excluded, including its 47 complete printed rows
and all damaged tail bytes. Every one of the 98,913 unresolved queries must
have a completed admitted full-database replay; absent output alone cannot
establish completed no-hit execution.

## Evidence Reviewed

1. The [archived source review](LEGACY_BLAST_RECOVERY_SOURCE_20260923.md)
   identifies the non-concatenated old-engine blastp path: one query is read,
   the search call returns, then tabular results are printed before the next
   query iteration. Re-read source lines 1900-1940, 2099 and 2213-2226 during
   this review. This explains why a subsequent query block is evidence for
   completion of printing the preceding one; it does not validate the last
   observed query or unobserved queries.
2. Rehashed the archived source, release archive, published MD5 file and
   installed binary. All SHA256 values agree with the earlier source review.
   The archived and installed 4,591,820-byte `blastall` files are byte-identical,
   SHA256 `34b91fedd2e7858a3478e83758e3f5339c09f96f543d9338488ea4b84628fb75`.
   No downloaded program was installed or executed. Source provenance does
   not establish a reproducible compiler/linker history.
3. The original exhaustive row audit, SHA256
   `a0d024ff33544b64eafcb7ebe88ba3948ad0cffe8aabcb1e4ef2ce300cc71bc9`,
   validates numerical/alignment structure and increasing query blocks.
   It records zero cutoff violations and zero logged-failed queries with
   outgoing hits. Its 44 failures are not treated as successful no-hit searches.
4. [Diagnostic 22055](qfo_blast_replay_comparison_22055.json), SHA256
   `9d49960fd22aaa1b538dcae728b06d0a67802d8d3d1f4296b31c459a98f702e1`,
   reproduces three complete earlier blocks exactly (20, 20 and 126 HSPs),
   the known failure with unchanged diagnostics, and a superset of the
   incomplete boundary block (73 versus 47 rows). Combined and isolated
   runs agree. This is limited sampling, not exhaustive native replay.
5. Read-only recheck **22148 completed 0:0 in 1:25**, using frozen executor
   `d61769462e67151482f2e817dc9926b26cbb33b1`. Its report
   [retained report](qfo_blast_prefix_recheck_22148.json) (originally
   `benchmarks/work/qfo_blast_prefix_recheck_22148.json`) has SHA256
   `315913c456280b55a5c591b7e4c138b5af58fc3c62a783fadde056473362d0d6`.
   Every indexed block hash/newline count matches; ranges remain contiguous.
   The complete 33,141,805,056-byte file, including the 966-byte damaged tail,
   still has SHA256
   `41063f82a09dd6b8ef4b4bbec783d0643e09f44279011dfb2687a1d43fe76f6d`.
   Descriptor/path identity remained unchanged during the scan.
6. The [frozen exhaustive partition](qfo_blast_recovery_partition_20260923.json)
   accounts for all 984,137 original queries: 885,224 proposed retained,
   98,912 absent and one incomplete boundary query. The replay set is not a
   smaller search database; all batches search the same full database.

## Conditions Still Required

Before production copying, bind this review and all listed evidence by
hash; require the completed recheck's scheduler/executor identity and fresh
file checks. Require the complete 20-batch panel's independent admission,
query-order/disjoint-coverage checks, exact database/runtime identity and
explicit failure/no-hit dispositions. Revalidate retained block hashes
while copying in original query order; never concatenate prefix then replay.

After copying, independently run the full numerical/alignment/diagnostic
audit and exhaustive query-disposition accounting. Only that separate
whole-search admission may reconnect downstream conversion and clustering.
Record interrupted and incremental costs separately as uncontrolled recovery
costs, not successful or matched comparative timing.

If any condition fails, preserve the failed artifact and investigate; this
review does not permit weakening the checks, automatic retries or silent
omissions. Original partial files and held job 21746 remain untouched.

## Residual Limitation

Historical disk-write ordering and unobserved corruption cannot be established
from a post-interruption file alone. The decision relies on ordinary sequential
output semantics, extensive structural and byte checks, official binary identity
and bounded diagnostic agreement. That uncertainty is retained explicitly.
An assertion of bit-for-bit equivalence to an uninterrupted full run would
require stronger evidence and is not made.
