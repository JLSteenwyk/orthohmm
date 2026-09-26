# Post-Failure Memcheck Artifact Integrity

Read-only checks of diagnostic 22164 completed using system Python with site
initialization disabled. No NumPy, igraph, Leiden, HMM refinement, optimizer or
runtime probe was executed. The [report](qfo_cpm_memcheck_artifacts_22164.json)
preserves all inspected file records and the original return code 97.

The status file matches its previously recorded SHA256. Its 990 input record
entries, nested runtime file records, child log and XML resolve to 1,053
distinct files, all matching their recorded identities. After checking the
output and metadata, 1,056 distinct files pass the final identity check.
Duplicate references are deduplicated only for hashing; conflicting identities
would be rejected.

The saved refinement metadata matches the successful native refinement except
for its output path. The output partition has identical bytes and SHA256
`f4c6f1973bc9636828baf1e6d9be3f416a18fc8ada5302fb502c082495fc1811`.
Independent streaming membership validation confirms all 984,137 genes appear
exactly once in 390,845 groups, with no foreign or duplicate members.

```bash
/usr/bin/python3 -S -m benchmark_tools.audit_cpm_memcheck_artifacts \
  --root . --output /tmp/qfo-cpm-memcheck-artifacts-new.json
python -m pytest -q tests/unit/test_audit_cpm_memcheck_artifacts.py \
  tests/unit/test_audit_failed_recovery_refinement.py
```

Six new artifact-audit tests plus ten existing partition/readback tests pass.
This completes the saved-file and partition checks left unfinished when the
diagnostic runner stopped at nonzero exit. It does not execute the runtime
probe omitted by that runner, establish memory safety, identify the SIGSEGV
cause, or admit high-CPM accuracy. The failed diagnostic and scientific
admission remain failed. No dependency was released or new inference launched.

The preceding startup-control prose was also corrected: its cwd and PYTHONPATH
were the checkpoint-recovery directory, not the original replay directory.
The raw control report already records that actual path. Its no-import result
remains useful but must not be described as an exact environment match.
