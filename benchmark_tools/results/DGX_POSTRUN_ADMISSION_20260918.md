# DGX Post-Run Admission

## Current Evidence

Final-task update: the controller subsequently reported21656_26 COMPLETED,
exit0:0, elapsed01:18:31,20CPUs/96GiB on spark-7ff0 with zero restarts.
Two sacct retries failed because localhost:6819 refused the accounting
connection. A later controller query no longer found the handle. The
successful response fields and later raw failure are retained separately in
dgx_final_task_controller_20260918.json; this is not a fresh all27 sacct
snapshot. Do not infer scheduler timezone from its unqualified timestamps.
Combined with the earlier completed0-25 snapshot, this establishes terminal
scheduler status for the timing sequence, not scientific resource admission.
Post-run evidence transfer and validation can now proceed; no rerun started.

`dgx_scheduler_progress_20260918.json` retains a timestamped `sacct`
snapshot for all 27 tasks in array21656. Tasks0-25 are COMPLETED with
exit0:0; task26 remains RUNNING. An independent parse verified unique
indices0-26,20 allocated CPUs,96GiB requested memory and spark-7ff0 for
every row, and no overlap between adjacent scheduler start/end intervals.
This is scheduler evidence, not admission of native output or resource data.

Every row reports TotalCPU00:00:00. Treat this accounting field as unusable
for CPU-efficiency calculations, not zero consumption. Do not substitute
allocated CPUs multiplied by elapsed time as observed CPU use. Scheduler
elapsed includes preparation and verification outside the native timer.

The earlier local run00 host replay reproduces its recorded inconclusive
status:22 unmatched identity events, all kworker-named, across12 intervals.
Names do not authenticate kernel threads. The maximum observed persistent
foreign load of0.007534 cores does not measure disappeared or unsampled
processes and cannot establish a controlled workload. This diagnostic
does not establish the status of the other26 runs.

## Required Validation

Post-run progress:207 files (815,526,081bytes) of metadata and measurement
evidence were transferred and checked with a checksum-mode rsync dry-run.
dgx_completed_panel_review_20260918.json records successful metadata,
payload-hash and host-replay checks for all27 runs. All27 host summaries
remain inconclusive:1,908 unmatched identity events across1,078 inconclusive
intervals. The other1,055 intervals report no large persistent competitor,
not proof of no competition.1,710 unmatched events are kworker-named, not
authenticated kernel threads. Native outputs and resource-accounting replay
remain unvalidated; no scientific timings have been admitted.

After the final task is authoritatively terminal:

1. Preserve all27 original outputs, metadata, raw host/resource samples and
   accounting. Transfer evidence only after the timing sequence finishes;
   retain source paths, sizes and content hashes.
2. Run the existing frozen metadata-contract audit on all27 indices. It
   checks command/input/runtime assertions and scheduler identity but is
   explicitly not a complete scientific timing admission.
3. Validate native output completeness and membership against frozen input
   inventories using validate_scaling_outputs.py. Keep native inference
   timing separate from preparation, verification, conversion and scoring.
4. Replay resource samples, verify monotonic clock domains and command
   bracketing, and reconcile collector accounting with the GNU-time
   companion. Report each memory statistic with its actual scope; do not
   interchange sampled aggregate RSS, cgroup memory and per-process maxRSS.
5. Review host intervals across all27 runs, including sampling failures,
   unmatched identities and actual foreign CPU use. Preserve inconclusive
   classifications. Any future enhanced observer must be versioned and
   calibrated separately; it cannot retroactively supply missing fields.
6. Determine reuse or selective reruns from the complete evidence and
   frozen criteria. Preserve failed/inconclusive runs and document any
   rerun decision. Do not silently label the whole panel controlled merely
   because Slurm allocations were exclusive or sequential.

No timing ranking, controlled-workload claim, or scientific resource table
is admitted by this document. Existing expensive results remain preserved;
no reruns or changes to the deployed collector have been initiated.
