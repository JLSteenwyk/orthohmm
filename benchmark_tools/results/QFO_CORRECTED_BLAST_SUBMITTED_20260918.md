# Corrected OrthoMCL Search Submission

Submitted job **21713** from detached executor
`benchmarks/work/publication_qfo_corrected_blast_v1`, revision
`b0c7128d464e2abe5bf1fd3646c0378625bbb0de`.
Scheduler inspection reports **PENDING, Resources**, with 180 CPUs,
900 GiB, 14-day limit on bizon, no requeue and zero restarts. Search has
not started and no start-time estimate is available. No competing job was
stopped or modified.

The standalone real preflight and installed-engine smoke passed before
submission; 19 focused tests and batch syntax validation passed. The
scheduled job must recheck all frozen records before native execution.
It binds prepared manifest
`3106012dc42c053d42e7f9d8d08532168d5826aa6812bab4b4c232f900e3a8ff`
and runtime inventory
`9ce46b5329f34384b03980b62fbfe9522f244dc227ce6506071ce7000b154e42`.

Monitor `benchmarks/work/qfo_corrected_blast_21713.log` and, after startup,
`benchmarks/results/qfo_corrected_orthomcl_v1/search_execution/status.json`.
This job executes formatdb and legacy BLAST only. Sequence-specific errors,
database completeness and query coverage require explicit auditing even
if both commands exit zero. BPO/index construction, native OrthoMCL
inference, final-group conversion and QfO scoring remain required and
are not automatically admitted or launched by this job.

These are shared-host stage measurements, not matched timing evidence.
Existing corrected comparator jobs, original-release factorial assessment
21711/admission 21712 and dedicated DGX timing remain separate.
