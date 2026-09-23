# Recovery Merge Component

`merge_blast_recovery_blocks.py` implements the streaming part of the frozen
[recovery protocol](QFO_BLAST_RECOVERY_PROTOCOL_20260923.md). It has no
production CLI and has not read, copied or merged the real interrupted BLAST
table. No scheduler dependency or held downstream job was changed.

## Implemented Checks

- Traverse the exhaustive original query order, not prefix-then-replay order.
- Require exactly one retained block for every non-replay query, contiguous
  prefix byte ranges and the exact authorized cutoff; reject the old final
  incomplete block.
- Preserve replay query order and permit an absent output block only for a
  query in the replay universe. Completion and no-hit/failure classification
  must already be proven by the independent complete-panel admission.
- Reject duplicate original/replay IDs, repeated output query blocks,
  unconsumed inventories and unknown source handles.
- Copy bounded blocks, checking every row's query ID, newline, absence of NUL,
  twelve fields, row count and exact block SHA256. Numeric/alignment validation
  is deliberately still required on the completed candidate.
- Check source inode/device/size/timestamps before and after copying, preserve
  the source files, fsync candidate bytes before rename, then fsync the output
  directory. A renamed candidate is not an admission or execution-success marker.

Twenty synthetic tests pass, including interleaved missing queries, boundary
replacement, excluded broken tails, wrong hashes/counts/IDs, malformed rows,
late iterator failure, changed sources, duplicate output queries and fsync
failure. Tests compare the exact expected merged bytes and digest.

```bash
python -m pytest -q tests/test_merge_blast_recovery_blocks.py
```

## Required Before Production

The driver still needs to bind this component to the hash-checked original
query partition, per-query inventories, complete twenty-batch admission and
an explicit prefix-reuse decision. It must verify full input hashes and each
retained range again, retain diagnostic/failure provenance, record execution
costs and write a durable terminal report. The component itself is not that
driver and cannot establish missing-query completion or prefix durability.

After merge, independently audit the entire table and exhaustive query
dispositions before BPO conversion or any scientific admission. No automatic
reuse, native retry, source-to-binary equivalence claim or final OrthoMCL score
is authorized by these tests. Original interrupted artifacts remain unchanged.
