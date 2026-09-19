# Corrected SwissTrees Strata Submission

Job21896was submitted held at2026-09-19 01:17:27EDT, then released after
checking its exact executor command,2CPU/64GiB allocation request onbizon,
4-hour limit, no-requeue setting and dependencies:
`afterok:21894:21736_0`. It remains pending, not a completed result.
21894produces the complete corrected factorial admission inventory;
21736_0independently admits the full OrthoFinder3.1.5 assessment. The
sequence-only sibling21736_1is not substituted.

## Frozen Execution

Committed and pushed executor revision:
`ba75fdcf27b8f58a425d674470ba8d889c81415a`, retained as a clean detached
worktree at `benchmarks/work/qfo_corrected_strata_executor_ba75fdc`.
Its HEAD and tracked benchmark-tool source cleanliness were checked before
release and are checked again by the batch.

- [Batch source](qfo_corrected_swiss_strata_batch_20260919.sh):
  `83fc3971ee41e5871b0b623594a9c260a84e6d421410d88c2cedf8694b714cb7`.
- Unchanged `run_corrected_swiss_strata.py`:
  `a4cc686b9c274137ad78b5814e493b2651c5f6d56c2ddbf4f9ae627532eedac0`.
- [Frozen protocol](CORRECTED_SWISS_SEQUENCE_STRATA_PROTOCOL_20260918.md):
  `2c87ba1df8ee39dfc325eefcdf374e713ba8fab1279a0330c3ea600ab637da54`.
- Frozen input-only descriptor inventory:
  `912363276fa5456f11aeff9f398fce71aec8d377c8c8eda7b144105945e554f1`.

All17driver dependency hashes and the descriptor/protocol hashes were
rechecked before submission. Python3.10.13 and NumPy2.2.6 are explicitly
required and recorded to the job log at execution. Hash seed0 and BLAS/OpenMP
thread limits1 are set. This local analysis is not a timing benchmark and
does not access the DGX.

The batch reads the final admitted files only after both scheduler
prerequisites succeed, records their actual hashes and passes them into the
existing strict driver. The complete factorial inventory comes from
`benchmarks/work/qfo_corrected_factorial_uncertainty_v1/scores/manifest.json`.
Full OrthoFinder admission comes from
`benchmarks/work/qfo_corrected_orthofinder_full_assessment_admission_20260918.json`.
The driver reconstructs raw family counts, checks provenance and exact
membership, selects p1c0r0/p1c1r1/full OrthoFinder, recomputes the frozen bins,
and applies100000paired replicates with seed20260924and27primary endpoints.
No partial input, historical prediction, changed method or replacement
stratum is accepted. Current low/high entropy bins are not new annotations
of protein fragments, domain architecture or evolutionary divergence.

## Verification and Remaining Work

`bash -n` passed. The driver, bootstrap, descriptor, corrected factorial
count and corrected comparator count test suites passed73tests in0.84seconds.
An initial test command used a nonexistent test filename and ran no tests;
the corrected invocation included `test_corrected_swiss_sequence_strata.py`
and completed successfully. No full-data stratified outcome was evaluated.

Fresh result: `benchmarks/work/qfo_corrected_swiss_primary_strata_v1.json`.
Log: `benchmarks/work/qfo_corrected_swiss_strata_21896.log`.
Missing/symlinked prerequisites, changed source/runtime and occupied output
paths fail rather than being bypassed. Failed scheduler prerequisites keep
this job from starting. Terminal validation, independent result checks,
figures/tables and manuscript integration remain required. All-method and
secondary-stratum displays are separate unfinished work; this submission
does not establish accuracy advantages, causal mechanisms or readiness.
