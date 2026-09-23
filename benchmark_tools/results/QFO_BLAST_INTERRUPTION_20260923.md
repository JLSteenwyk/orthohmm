# Corrected QfO BLAST Interruption

## Evidence

On September 23, the controller still reported job 21713 RUNNING, but this
was not a live computation. Direct host inspection established:

- `uptime -s`: 2026-09-23 09:19:58 EDT; Slurm node BootTime was
  2026-09-23T09:19:57 and SlurmdStartTime was 09:20:31.
- No `blastall` process or job-specific worker was present. The only
  observed slurmstepd was the service's `infinity` process.
- `scontrol listpids 21713` returned exit 1, "JobId=21713 does not exist
  on node (null)."
- The slurmd startup journal recorded removal of the vestigial
  `/var/spool/slurm/slurmd/job21713/slurm_script` at 09:20:31.
- `work/all.blast.partial` retained 33,141,805,056 bytes, last modified
  2026-09-22 15:58:23.012870347 EDT. This timestamp predates the reboot;
  it does not establish the exact interruption time or its cause.
- The retained execution status says `running_blast` but contains only
  a BLAST start timestamp, not a finish or exit code. `blast.time.txt`
  and the batch log are both empty. These stale records are preserved,
  not rewritten as successful execution.

Paths above are relative to
`benchmarks/results/qfo_corrected_orthomcl_v1`, except the explicitly
absolute Slurm path. Observations are transcribed from command output;
they are not a complete archived host journal.

SHA-256 fingerprints taken after confirming the worker absent:

| Retained file | SHA-256 |
| --- | --- |
| `work/all.blast.partial` | `41063f82a09dd6b8ef4b4bbec783d0643e09f44279011dfb2687a1d43fe76f6d` |
| `search_execution/status.json` | `1e395daea118a44d482af9107d059085c655382c74ea07bbbbc85d7f622a460a` |
| `search_execution/blast.log` | `876f55697e53bf937ef3a0b4f5f021ab7350853752251fee814f7002ef8b0617` |

## Recovery Actions

Held downstream admission job 21746 with `scontrol hold 21746`, then
cleared the nonexistent BLAST allocation with `scancel 21713`. Accounting
now records `CANCELLED by 1000`, ending 2026-09-23T09:39:22, with controller
elapsed time 4-01:29:59. The displayed 0:0 is NOT a native BLAST success;
that elapsed time is NOT verified computation time after the interruption.

The remaining OrthoMCL dependencies were left intact behind held job
21746. No partial table was renamed, truncated, converted, or admitted.
The queued corrected FastOMA job 21740 then started; its wrapper,
Nextflow, and containerized input-check process were directly observed.
Unrelated processes and services were left alone.

## Outstanding Work

The corrected OrthoMCL result is incomplete. The frozen runner deliberately
rejects existing outputs and does not implement resumable BLAST. Do not
relaunch it against this directory or infer completed queries solely from
the last output row: no-hit queries and interrupted query output need
explicit handling. Recovery must preserve these files, validate any
reusable query results against full-database semantics, retain failure
diagnostics, and record a separate execution plan before further work.
Otherwise a fresh, separately recorded search is required. Keep 21746
held until a reviewed recovery replaces or reconnects that chain.

The scientific configuration and endpoints are unchanged. No final
OrthoMCL score or comparable timing is claimed from this interrupted run.

## Read-only Byte Audit

The subsequent full-file scan using `audit_interrupted_blast_bytes.py`
reproduced the preserved SHA-256 and found:

- 366,012,663 newline bytes; the last is at zero-based offset
  33,141,804,089. Newline count is not a validated HSP-row count.
- An 11-byte fragment follows that newline, then exactly 955 NUL bytes
  at offsets 33,141,804,101 through 33,141,805,055 inclusive.
- All NUL bytes form one trailing run. None were found in earlier bytes.
- The file does not end in a newline. A bounded tail inspection showed
  the truncated final row starts `tr|F1MRE1|F`; preceding rows are for
  query `tr|F1MRE1|F1MRE1_BOVIN`.

[Machine-readable inventory](qfo_blast_interrupted_bytes_20260923.json)
explicitly sets search admission and reuse authorization false. The audit
checks file identity/size/timestamps before and after scanning and never
edits the input. Its 33 focused tests passed, including chunk boundaries,
interior versus trailing NULs, empty inputs, no-overwrite behavior, and
input mutation rejection.

This narrows byte corruption to the observed tail but does not certify
earlier rows, query completion, or durable write ordering. No rows were
salvaged. The installed legacy documentation does not establish the needed
query-completion guarantees. Attempts to retrieve source listings from
`https://ftp.ncbi.nlm.nih.gov/blast/executables/legacy/` and its `2.2.13/`
child returned HTTP 404. Modern BLAST+ behavior was not substituted as
proof of legacy behavior. Next: establish a validated query-completion
boundary and replay strategy, or use a fresh separately recorded search.

## Retained Row and Query-Order Audit Completed

Job 22030 completed with exit 0:0 in 00:44:28 on September 23. The
read-only scan applied the existing numerical, alignment-accounting,
identifier, coordinate, and query/subject-block checks to every complete
row before the frozen damaged-tail boundary. It also confirmed strictly
increasing FASTA query ordinals between blocks and reproduced the original
full-file SHA-256, including the excluded tail.

| Observation | Count |
| --- | ---: |
| Structurally validated HSP rows | 366,012,663 |
| Distinct directed query/subject blocks | 202,743,309 |
| Observed query blocks (including incomplete final query) | 885,225 |
| Queries with self hits | 885,225 |
| Subjects with hits | 958,574 |
| Input queries without observed hit rows | 98,912 |
| Logged failed queries | 44 |
| Logged failed queries with outgoing hits | 0 |
| Logged failed queries with incoming hits | 16 |
| HSP rows exceeding the configured 1e-5 cutoff | 0 |

The final query `tr|F1MRE1|F1MRE1_BOVIN` begins at zero-based byte offset
33,141,800,005. Its 4,085 bytes of complete rows and the following 966
damaged-tail bytes are NOT a complete query result. All of that query must
be excluded from any proposed retained prefix. This leaves 885,224 earlier
observed blocks as candidates for further recovery validation, not accepted
results. A conservative replay set consisting of all absent queries plus
the final query would contain 98,913 unique input queries. Absence includes
unprocessed queries and potentially true no-hit/failed queries; do not call
98,912 a native search-failure count.

Artifacts remain outside Git under
`benchmarks/work/qfo_blast_prefix_audit_20260923/`:

- `status.json` SHA-256:
  `a0d024ff33544b64eafcb7ebe88ba3948ad0cffe8aabcb1e4ef2ce300cc71bc9`.
- `query_blocks.jsonl`: 229,006,947 bytes, SHA-256
  `fc0c4cb91088308b64674adf657ff3df63e234429110bfbd263856395d575ad6`.

The report explicitly leaves admission and reuse false. The existing
44 diagnostic failures are incomplete-run observations, not the final
benchmark failure inventory. Native replay 22029 remains pending; no
retained prefix was copied, merged, or admitted. FastOMA 21740 remains
live and has advanced to `hog_rest` tasks.
