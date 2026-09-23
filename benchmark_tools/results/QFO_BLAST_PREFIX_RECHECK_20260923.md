# Retained Prefix Byte Recheck

Before the prefix evidence review, reread every block in the pinned original
885,225-block inventory and compare its SHA256 and newline count. Require
contiguous ranges, increasing original query ordinals, one final boundary
block and the original total-file SHA256 including the damaged tail. Only
885,224 blocks ending at byte 33,141,800,005 are retention candidates; the
entire last query remains excluded, including its 47 complete printed rows.

The new reader uses bounded 8-MiB chunks rather than reparsing hundreds of
millions of numeric rows. Structural/alignment checks are supplied by the
hash-bound earlier exhaustive audit; this is an independent byte-range and
full-file identity recheck, not a second numerical validator. Recheck the
earlier report, FASTA, diagnostics, block inventory and validator/source
identities before and after. Check open-descriptor and path metadata against
replacement or mutation during the scan. Preserve the original file read-only.

Run once from a frozen executor with 2 CPUs, 8 GiB, a 4-hour limit and no
automatic requeue using `qfo_blast_prefix_recheck_20260923.sh`. Output is a
new job-specific report, never a modification of the earlier evidence.
Thirty scan/merge tests pass, including altered ranges, rows, hashes, order,
NULs, boundary, tail, truncation, duplicate blocks and failed fsync.

The result deliberately retains `reuse_authorized=false` and
`search_admitted=false`. Combine a successful recheck with the archived
sequential-query source path, official-release binary identity, diagnostic
agreement and exact query partition in the separate evidence review. It
cannot prove historical filesystem durability or reproducible source-to-binary
compilation. All twenty replay admissions and final merged-table validation
remain required; held job 21746 and downstream inference stay untouched.
