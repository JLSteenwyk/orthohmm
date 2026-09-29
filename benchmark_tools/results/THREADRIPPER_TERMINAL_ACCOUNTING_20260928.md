# Terminal Scheduler Accounting Gap

Fresh read-only `sacct` queries confirm jobs 22367-22369 and their batch/native
steps all COMPLETED 0:0. Allocation durations are 145, 146 and 146 seconds;
native-step durations are 12 seconds each. These are coarse scheduler elapsed
times for prior 16-gene fixtures, not comparable scientific runtime results.

All nine rows have blank MaxRSS/MaxVMSize and TotalCPU `00:00:00`. Current
`scontrol show config` reports `JobAcctGatherType = (null)` with
`AccountingStorageType = accounting_storage/slurmdbd`. Storage of terminal
job state does not imply collection of resource usage. The
[official accounting documentation](https://slurm.schedmd.com/accounting.html)
distinguishes accounting gathering from storage; the
[sacct documentation](https://slurm.schedmd.com/sacct.html) warns that missing
resource data can appear as zero. Current configuration is not proof of the
historical configuration, but these records themselves lack memory evidence.

The [captured commands and assessment](threadripper_terminal_accounting_20260928.json)
retain raw stdout/stderr, timestamps and inspector identity. The standalone
inspector treats all nine CPU and memory usage rows as unavailable under the
observed absent collector, never measured zero. Nine tests pass, including
missing/duplicate/wrong-job records and non-admission even when fields exist.
This inspector is diagnostic; it is not integrated production authorization.

## Consequences For Timing

- Native and job cgroup measurements already retained remain valid only for
  their recorded windows. This observation does not invalidate those reads.
- The reported job peak ends at report generation. No terminal Slurm value
  can extend it through subsequent validation and teardown for these fixtures.
- Do not sum step MaxRSS or relabel an in-job read as a complete-job peak.
- Before production, validate an external observation path that captures
  final job counters after owned work exits and before cgroup removal, or
  obtain and validate an appropriately configured accounting mechanism.
  Missing final observations must remain unresolved, not inferred from zeros.
- No scheduler configuration, daemon or unrelated workload was changed. No
  job was submitted. An administrator-mediated configuration change would
  require separate coordination and validation; it was not attempted here.

Whole-run environmental policy, native handoff validation, full-scale observer
overhead and a quiet window remain separate prerequisites. The original
27-run identities and Threadripper-only instruction are unchanged.

```sh
python -B -m benchmark_tools.inspect_terminal_accounting --job 22367 --job 22368 --job 22369 --output /tmp/threadripper-terminal-accounting-new.json
```

This polls retained jobs, not new timing runs. Use an unused output path.
