# OrthoMCL Search Table Audit

## Scope

`audit_orthomcl_search_table.py` adds a streaming, pre-BPO structural and query
coverage check. It does not modify the frozen search engine, conversion code,
masking defaults, thresholds, or queued corrected BLAST job 21713.

Checks include twelve-column protein m8 format; known query/subject identities;
finite, nonnegative scores and E-values; bounded identity and alignment counts;
coordinates within input sequences; consistency of alignment spans, gaps,
mismatches and two-decimal identity; and contiguous query/subject blocks.
Repeated adjacent HSPs are supported. Reopening a closed block is rejected
because the existing BPO converter combines adjacent HSPs into one pair entry.
Such a failure requires investigation, not automatic sorting or repair.

The audit reports query, subject and self-hit coverage, diagnostic categories,
logged failed queries with incoming/outgoing hits, and all input IDs lacking
reported query hits. No-hit queries without a logged failure are explicitly
separate. HSPs above the nominal E-value cutoff are counted, not filtered out.
File and helper hashes are checked before/after processing; existing output
reports cannot be overwritten.

## Native Smoke Evidence

The installed legacy BLAST 2.2.13 and formatdb were run with one thread in a
fresh `benchmarks/work/orthomcl_table_smoke_20260918` directory. The existing
opt-in native-engine test now checks both two positive-control sequences and
the same sequences plus a two-residue query. Raw inputs, database files, m8
tables and logs are retained there, separate from production inputs/results.

| Fixture | Inputs | Directed pairs | Queries with hits | Logged failed queries |
| --- | ---: | ---: | ---: | ---: |
| Positive controls | 2 | 4 | 2 | 0 |
| Positive controls plus short query | 3 | 4 | 2 | 1 |

Both native commands exited zero in both fixtures. The short query produced
one setup warning and one short-query error, correctly counted as one failed
query, with no incoming or outgoing hits. This directly demonstrates why a
zero BLAST process exit is insufficient evidence of per-query success.

Saved reports are `orthomcl_table_smoke_success_20260918.json` and
`orthomcl_table_smoke_failure_20260918.json`. These are tiny native format and
failure-handling tests, not a full corrected database or biological validation.

## Remaining Gates

Production use must first bind completed job 21713, the frozen search executor,
execution report, input parity and runtime manifests. The formatted database
must be independently checked against the complete corrected sequence set.
The full table can then be audited and its diagnostic query failures reviewed.
BPO/index validation, native final clustering, final-group conversion and
reference-impact/scoring checks remain separate unfinished steps. Historical
failed-query counts must not be transferred to the corrected release.

No search, accuracy, timing or publication admission is asserted by the new
table auditor. Existing native diagnostics fail closed on unrecognized output;
the full production log may require additional documented diagnostic handling.
