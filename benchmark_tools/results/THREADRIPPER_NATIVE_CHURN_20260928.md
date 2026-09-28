# Native Process Churn Is Not Foreign Competition

A review of the new whole-run process-policy path found an applicability bug:
it rejected identity changes for every persistent PID, including native children
observed inside the benchmark's own job in both snapshots. Three added tests
failed before the fix: native executable/name changes, later PID replacement
with reset counters, and movement between steps inside the same job. Normal
birth and exit were already excluded by the existing scoped CPU diagnostic.

For the typed v2 policy only, non-observer processes observed inside the job
at both endpoints now retain identity transitions in
`observed_in_job_identity_changes` without classifying them as external work.
A later creation timestamp identifies a different process, so its reset CPU
counter is not compared as if it belonged to the earlier process. A decreasing
creation timestamp fails. A counter decrease for the same process still fails.
Historical v1 policy semantics and the frozen shared CPU diagnostic are unchanged.

This is not a blanket PID exemption. The observer must retain its identity and
membership. Foreign-to-job and job-to-foreign transitions, unreviewed external
processes, missing type evidence and collection errors still fail. CPU use of
reviewed kernel threads remains visible. Native CPU and memory measurements
remain in their separate cgroup/lineage records; no resource amount is subtracted.

273 focused tests pass, including the previously failing cases, scope-boundary
and counter controls, and a complete process-stream test that accepts native
exec but rejects movement outside the job. These are synthetic regression tests,
not a native Slurm or full-scale observation validation. No existing timing was
readmitted, output rescored or numerical threshold relaxed.

Endpoint membership cannot establish that a process remained in the job between
observations. In particular, uncertain collection errors are not waived merely
because an earlier row belonged to the job. Short-lived native processes may
still cause incomplete observations; resolving such gaps requires better evidence,
not assuming they are harmless. Native/full-scale validation and the reviewed
quiet-window policy remain required before production timing.
