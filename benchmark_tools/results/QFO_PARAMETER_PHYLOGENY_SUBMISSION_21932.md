# Corrected QfO Candidate Phylogeny Submission

Submitted 2026-09-19 at 08:23 EDT as array **21932**, indices 0-3,
concurrency one, ordered norm_low, norm_high, margin_low, margin_high.
All four tasks were PENDING at the first scheduler check. No native result
or accuracy result is claimed.

## Frozen Execution

- Executor commit: `aa8c0e1937b898a9da83c69bf36ec342a4e04b89` (pushed).
- Detached executor: `benchmarks/work/publication_qfo_parameter_phylogeny_v1`.
- Runner: `benchmark_tools/run_qfo_parameter_phylogeny.py`, SHA-256
  `195b176a2b90e08f42066d56a7a108f177aba59eb5dd519a235715dbf3d886c5`.
- Batch: `benchmark_tools/results/qfo_parameter_phylogeny_batch_20260919.sh`,
  SHA-256 `0fea95574656865b56f2b7ba94c84f009b9a35c85f4707d8bc55f40882f70acc`.
- Candidate admission: retained 21929 JSON, SHA-256
  `77759ee0f733242d6d4c368f44c2c1fd84ef3b1ebbb679494acce7fab68f0141`.
- Baseline native admission: retained 21764 JSON, SHA-256
  `c8a289ac8128711da6c4e93654ce6953e48854b8d21b2adfea074d3feb8c70aa`.

The submission used the batch inside the detached executor, passing that
executor's absolute path and exact commit as positional arguments. Slurm
confirms 32 CPUs, 192G RAM, 24 hours per task, node bizon, no requeue and
array throttle one. The initial scheduler reason was
`Nodes_required_for_job_are_DOWN,_DRAINED_or_reserved_for_jobs_in_higher_priority_partitions`;
this reason does not establish that bizon itself is down. OrthoMCL 21713
was concurrently confirmed RUNNING with 180 CPUs. No resource limits or
unrelated jobs were changed to accelerate scheduling.

## Outputs and Remaining Gates

Outputs will be under `benchmarks/results/qfo_parameter_phylogeny_v1/<arm>`;
batch logs and GNU time files are `benchmarks/work/qfo_parameter_phylogeny_21932_<index>.*`.
Preflight records bind inputs, source, helpers, exact command and tools.
Baseline artifacts and runtime are rechecked after inference. Successful
execution is explicitly pending independent native-output validation.

Next: validate terminal accounting and native output inventories, admit and
convert native ortholog pairs, score the six frozen endpoints, and implement
the prespecified paired SwissTree contrasts. Do not infer family-level
TreeFam uncertainty from the pooled reference. The separate CPM variants
remain outstanding. Shared-host checkpoint-reusing runtime is incremental
execution evidence, not a matched end-to-end timing comparison.

Verification before deployment: 80 focused tests passed; `bash -n` passed.
These tests exercise admission selection, scheduler gates, execution/failure
handoff, working-directory restoration, refusal to overwrite existing
outputs and inherited command constraints. They are not an end-to-end test
of this new panel; scheduled execution and independent admission remain
necessary.
