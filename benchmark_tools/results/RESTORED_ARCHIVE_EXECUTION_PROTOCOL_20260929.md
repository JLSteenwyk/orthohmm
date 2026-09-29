# Restored Archive Full OrthoBench Execution

## Prospective Scope

Execute the restored local asset archive with the independently acquired public
OrthoBench input manifest and reconstructed base/installer Python. This tests
the archive-to-results path, a distinct reproduction requirement from completed
job 22376. Do not replace or rerun that job. This is same-host reproduction,
not independent biological validation, controlled timing or a public release.

Use `/tmp/orthohmm-restored-execution-assets-20260929-v2` and
`/tmp/orthohmm-public-input-toolkit-20260929/rebound/data.json`.
Preserve the failed first restoration and all prior receipts. Use the existing
integrated controller and job launcher without modifying their scientific
settings. Rebind the previous plan's file records to restored assets and newly
acquired input paths; verify content identities against the previous pins.
Retain exact base-runtime pins. Record the new plan digest before submission.

One fresh attempt, no checkpoint reuse or automatic retry; 32 CPUs, 128 GiB,
24-hour Slurm allocation on bizon. Do not stop other workloads or alter services.
Retain failure logs and incomplete outputs. Shared-host elapsed time is diagnostic
only. Do not submit if a matching attempt already exists or is running.

## Required Admission

Require scheduler success and all eight controller stages to succeed. Verify
the plan and dependencies after execution, freshly installed package payloads,
complete native outputs, independent readback and all 70 RefOG score objects.
Compare the full label-invariant root partition against admitted job 22376;
compare the four retained native pair/confidence/reconciliation/hierarchy TSVs
byte-for-byte. Report any differences, rather than relaxing acceptance or
changing the algorithm. Bind independent admission to this specific new job and
plan; do not change the job-specific validators for earlier completed runs.

The archive excludes the base interpreter, installer, OS libraries and raw
datasets. Successful execution does not establish hermeticity, cross-host
portability, redistribution rights or publication readiness.
