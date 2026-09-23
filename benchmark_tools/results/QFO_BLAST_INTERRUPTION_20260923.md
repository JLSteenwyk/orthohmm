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
