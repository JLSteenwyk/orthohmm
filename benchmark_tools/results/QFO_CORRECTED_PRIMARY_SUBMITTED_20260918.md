# Corrected QfO Primary Submission

## Execution

Submitted array **21706**, tasks 0-1 with concurrency one, from detached
executor `benchmarks/work/publication_qfo_corrected_primary_v1` at
`9d8b608ab9b4a152d67ceacc04fe49b4d7595788`.
The batch script is the committed
`benchmark_tools/results/qfo_corrected_primary_batch_20260918.sh`.
The command manifest is `qfo_corrected_primary_commands_20260918.json`,
SHA-256 `dbebd3a6915fddeb2b89e0e591a2ac5c798ee6bccec9925fef485260f41baa5a`.
The earlier preparation manifest remains immutable; this ledger records
authorization and submission after the pinned-launcher gate passed.

Scheduler inspection confirmed task 0 running on bizon with 32 CPUs,
192 GiB memory, 72-hour limit, no requeue and no restart. Its numeric
job ID is 21707. Task 1 is pending the array concurrency limit.
Task 0's native log confirms 78 FASTAs, high-sensitivity HMM search,
BLOSUM62, E-value 1e-4, Leiden CPM 0.1, and 32 CPUs; all-to-all search
has begun. The execution status records the exact manifest hash, native
arguments, frozen core working directory and pinned runner source.

Task 0 runs OrthoHMM high-sensitivity; task 1 runs full OrthoFinder 3.1.5.
Outputs are isolated under `benchmarks/results/qfo_corrected_primary_v1`.
Neither output is yet admitted or scored. Resources measured here will
be descriptive shared-host measurements, not dedicated timing results.
The separate DGX array is unchanged.

## Scope Still Open

Complete native validation, conversion and scoring for these runs; freeze
and execute the six other comparator rows and all eight corrected-input
factorial cells. Original-release results remain separately labeled.
No old-release search checkpoint or output group was reused for this run.
