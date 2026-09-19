# Corrected QfO Candidate Neighborhood Submission

Job 21927 was submitted September 19, 2026 and confirmed RUNNING on bizon
at 07:57:32 EDT, 14 seconds after start. It requests 2 CPUs, 64 GiB and four
hours, with no requeue. This is shared-host candidate preparation, not a
timing benchmark or completed accuracy evaluation.

The detached executor is
`benchmarks/work/publication_qfo_parameter_candidates_v1`, revision
`6405da57618cc703fc1ca0f51e450d3cf6e2a05f`. Its tracked analysis/scientific
files were clean before submission. The batch verifies that exact revision
and tracked source state before invoking Python. The submitted script is
the committed `benchmark_tools/results/qfo_parameter_candidates_batch_20260919.sh`
inside that executor. Submission cleared `LD_LIBRARY_PATH`, `LD_PRELOAD`,
`LD_AUDIT` and `PYTHONPATH`, and set `TMPDIR=/tmp`; the batch fixes hash seed
and numerical thread settings.

- Preparation source SHA-256:
  `e6e7d889fe0e0cc545ac9c01c84ae146ca80828324bb49e4eb9023c4e3263aeb`.
- Batch source SHA-256:
  `5aa182a5aee66c2a77162fc63b5d54f6ae6e98fe0958e1b68f835fbd5a0ed777`.
- Outputs: `benchmarks/results/qfo_parameter_candidates_v1/`.
- Log: `benchmarks/work/qfo_parameter_candidates_21927.log`.
- Incremental GNU time record:
  `benchmarks/work/qfo_parameter_candidates_21927.time`.

The driver revalidates the pinned corrected baseline and scientific runtime,
prepares the unchanged control followed by four threshold variants, and
requires exact control candidate/merge bytes before any variant. Candidate
content, trace consistency, applied parameters and partial failures are
retained. 69 focused tests and batch syntax checks passed before submission.
Full independent preparation admission, both CPM variants, native phylogeny,
pair conversion and scoring remain outstanding.

## Retained Submission Failure

Job 21926 FAILED with exit 1:0 after one second because the invocation passed
the incorrect revision string `6405da53bfbd89374c06fe08587d8b6fef020ae0`.
The batch's first revision comparison rejected it before Python/GNU time.
The proposed output directory did not exist after failure; no candidates or
scientific results were generated. The replacement corrected only the
invocation argument to the verified worktree hash above. No source/settings
change, outcome-based retry or deletion of failure evidence occurred.
Its scheduler accounting and `qfo_parameter_candidates_21926.log` remain
retained locally. This failed submission must remain in the execution history.
