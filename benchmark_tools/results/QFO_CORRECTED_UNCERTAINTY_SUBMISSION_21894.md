# Corrected Factorial Uncertainty Submission

Submitted21894on19September2026 at01:09:20EDT, initially held. Verified its
requested resources, exact executor command and all three dependencies, then
released it. It remains pending on successful completion of independent
score admissions21780,21784and21788. Earlier five cell receipts are already
available; all eight are required at execution, not just the new three.

## Frozen Execution

- Committed and pushed batch revision: `9dbf842a633258c8896bd4feb65f0494a4d6c9ff`.
- Detached executor: `benchmarks/work/qfo_corrected_uncertainty_executor_9dbf842`.
- [Batch source](qfo_corrected_factorial_uncertainty_batch_20260919.sh) SHA-256:
  `5acdabeba75c9476196bf908a04de46ec0a764cb47a0b05dc73faa1507611127`.
- Local hostbizon, partitiongpu,2CPUs,64GiB,4hours, no requeue.
- `afterok:21780:21784:21788`; no automatic bypass of failed prerequisites.
- Python `/home/bizon/anaconda3/bin/python`; hash seed0 and BLAS/OpenMP
  thread limits1. This is uncertainty computation, not comparative timing.

Executor HEAD and clean tracked benchmark-tool sources were checked before
release; the batch repeats those checks. Required protocol/implementation
hashes match the existing bootstrap's frozen constants:

| Artifact | SHA-256 |
| --- | --- |
| Corrected release protocol | `b3603bc9b51ce1708f49a0ee816eb02d3cb52c66fb00f07eeb0d66f840a03e1e` |
| Factorial protocol | `f8946e12cefcf84abbee0fb9492f240c05508e045efe00a3006304d34c1fd115` |
| Corrected count auditor | `375f2b4fd511c04e562827aebacb3b7522ca979f76b665652dd209940542f912` |
| Shared numerical engine | `be09876ae4de31b923818385bd3d40d8d5e7df177215fc45a2101c7d30d1c1fa` |

The batch captures hashes of the eight final admission receipts at execution,
then the existing exporter checks each receipt and its bound conversion.
The count auditor validates the complete frozen cell order, independently
admitted execution/raw-file bindings and reconstructed SwissTrees aggregate
against each admitted endpoint. Historical reference evidence anchors labels
only; historical prediction counts are not substituted. All checked inputs
and their hashes are carried into generated provenance.

The unchanged analysis uses100000paired family replicates, seed20260922,
18SwissTrees families and42multiplicity-adjusted endpoints, recomputing each
aggregate from resampled counts. No partial-factorial significance test,
parameter retuning or new endpoint is introduced.

## Outputs and Validation

Fresh output root: `benchmarks/work/qfo_corrected_factorial_uncertainty_v1`.
Expected outputs are `scores/manifest.json`, `scores/scores.tsv`,
`scores/scores.md`, `swiss_counts.json`, `swiss_bootstrap.json` and
`swiss_bootstrap.md`. Existing output roots are rejected. Missing/symlinked
admission files are rejected before export. Failures after partial work leave
that evidence in place; the batch does not overwrite it on retry.

Log: `benchmarks/work/qfo_corrected_factorial_uncertainty_21894.log`.
`bash -n` passed; the existing exporter, corrected count auditor and corrected
bootstrap unit suites passed45tests in0.41seconds. This does not validate a
full-data result: execution, count/result inspection, independent arithmetic
checks, manuscript integration and retained final provenance remain required.
No uncertainty outcome or publication readiness is claimed at submission.
