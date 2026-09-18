# Corrected Candidate Preparation And Admission

## Frozen Work

The candidate engine and independent auditor were already implemented and
frozen. This milestone adds dependency-aware batch handoffs, not a new
algorithm or a parameter search.

- Producer: `0a028a743fca7626b86376475f4c7fd438093717` at
  `benchmarks/work/publication_qfo_corrected_candidates_v1`.
- Auditor: `9082ccee291176b8883884d80e80ff4817053b86` at
  `benchmarks/work/publication_qfo_corrected_candidate_admission_v1`.

Both source worktrees were checked clean. The producer requires successful
replay admission21757 and its frozen source09dec901. It builds all four
`p0_c0`, `p0_c1`, `p1_c0`, `p1_c1` arms from admitted corrected inputs, with
the two profile conditions using the refined multipass and refined profile
partitions respectively. Profile-off retains the initial HMM search; it is
not a fully HMM-free comparator. Candidate-on uses frozen satellite_v2
settings and membership constraints. The producer creates the eight-cell
plan but does not execute reconciliation or score any cell.

The independent auditor checks terminal preparation accounting, source and
helper inventories, corrected984,137-protein/78-species scope, numeric
checkpoint and species ownership, all four candidate arms, membership
constraints, the eight-cell plan and the fixed parameter values. Recorded
candidate search-support values are not independently recomputed by this
audit. Native/replay partition inequality is retained, not rejected or
silently treated as equality.

## Batch Handoffs

`qfo_corrected_candidates_prepare_batch_20260918.sh` must follow
`afterok:21757`; `qfo_corrected_candidates_admit_batch_20260918.sh` follows
`afterany` of preparation. Each requests2CPUs/64GiB/24h onbizon with no
requeue. The extended wall allowance covers full-input preparation and
verification; it does not change the frozen algorithm or turn shared-host
incremental time into matched end-to-end efficiency evidence.

Both scripts verify exact40-character executor revisions, clean source
trees and numeric job IDs. They compute the relevant input-manifest hash
after dependency satisfaction and invoke only the frozen producer/auditor.
Fixed Python hash seed and single-threaded numerical libraries are retained.

Preparation output:
`benchmarks/results/qfo_corrected_factorial_v1/manifest.json`.
Admission output:
`benchmarks/work/qfo_corrected_candidate_admission_20260918.json`.
Neither output authorizes an accuracy claim; reconciliation, pair conversion
and independent QfO scoring/admission remain required.

## Verification

16 new batch tests exercise syntax, resources, exact handoff arguments,
runtime-computed checksums, numerical environment, and rejection of missing
input/arguments, wrong revision, dirty source and invalid job ID. Tests use
isolated temporary git repositories and stub child programs, never forged
production manifests. Combined candidate preparation/admission/content and
batch suite:80 passed in0.81s.

An actual frozen-producer preflight, with explicitly named preflight
allocation variables, rejected pending21757 before reading partial replay
evidence or creating candidate output. This was not a successful scheduled
candidate run. Submission identities are recorded below after confirmation.
