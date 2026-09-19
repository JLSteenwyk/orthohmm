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

## Independent Admission Array 21935

Submitted 2026-09-19 at 08:29 EDT from detached executor
`benchmarks/work/publication_qfo_parameter_native_admission_v1`, commit
`72eb401ab056772f7f5b85585c8f30dc33bc3c78` (pushed).
The exact executor path and commit were passed to the frozen batch.

- Admission source SHA-256:
  `9c3ffb3d606be0e2ce11ff1883bf05ce567e3dad789071382d6976fd172dc66f`.
- Batch `qfo_parameter_phylogeny_admit_batch_20260919.sh` SHA-256:
  `e9e7adcf3e6bf25d7b10bb7ac1f1a68441853ca60306de2dece282f58b6ef365`.
- Four tasks, concurrency one, 2 CPUs, 64G, 4 hours, bizon, no requeue.
- Dependency: `aftercorr:21932`; each admission task requires successful
  completion of the matching inference-array index.
- Reports: `benchmarks/work/qfo_parameter_native_admission_21935_<index>.json`.
- Logs: `benchmarks/work/qfo_parameter_native_admit_21935_<index>.log`.

Confirmed all four admission tasks PENDING (Dependency); task 0 raw job
identity is 21936, distinct from array identity 21935_0. The admission tool
checks both array selection and recorded raw execution identity. Focused
admission/runner/native-validator suite: 114 passed in 8.95 seconds; batch
shell syntax passed. No report or native-output success is claimed before
these jobs execute. Mapping, scoring and paired uncertainty remain separate
downstream requirements.

## Native Pair Conversion Array 21939

Submitted 2026-09-19 at 08:37 EDT from detached executor
`benchmarks/work/publication_qfo_parameter_pairs_v1`, commit
`0401860cede4b1b73e5f6f9d459e0a0daabb1c6f` (pushed), passing the exact
executor path and commit to its frozen batch.

- Converter SHA-256:
  `ef209f6a68e70312cb7cfdfcac17c2f67b02bc2fb5744e4f457b8ad8b57e61b9`.
- Batch `qfo_parameter_pairs_batch_20260919.sh` SHA-256:
  `2a92a5891a58ab8769d5cdf5dfaa04a767a74bf145c66e7ef79296af7b8100c8`.
- Four tasks, concurrency one, 2 CPUs, 64G, 4 hours, bizon, no requeue.
- Corresponding-task dependency: `aftercorr:21935`.
- Outputs: `benchmarks/results/qfo_parameter_pairs_v1/<variant>/`.
- Logs: `benchmarks/work/qfo_parameter_pairs_21939_<index>.log`.

All four tasks confirmed PENDING (Dependency); task 0 has raw job ID 21940.
The converter reruns the pinned independent admission in a subprocess and
requires exact report equality before writing normalized native pairs. Any
unexpected reference-mapping loss fails conversion, preserving the partial
files and mapping counts. Tests exercise real conversion (not clique
expansion), mapping losses, mismatched fresh admission, output preservation
and refusal to retry existing output directories: 92 focused tests passed
in 0.54 seconds, and batch shell syntax passed.

No pair conversion has executed yet. Official assessment, score admission
and the frozen paired uncertainty analysis remain required downstream.

## Six-Endpoint Assessment Array 21944

Submitted 2026-09-19 at 08:43 EDT from detached executor
`benchmarks/work/publication_qfo_parameter_assessment_v1`, commit
`f5a42ab3fc6ee5076e9dfb2b239b63b131aa0d73` (pushed), passing the exact
executor path and commit to its frozen batch.

- Runner SHA-256:
  `e06867181061a19f3e2ef10054653e9f5798cd7b3be27bb3a3ca73da6972e60b`.
- Batch `qfo_parameter_assessment_batch_20260919.sh` SHA-256:
  `7bd70c2c3e8a83fae6b2ae6fb11e9b4165bfe8d1a5eadf9490d30145b222e161`.
- Four tasks, concurrency one, 8 CPUs, 96G, 24 hours, bizon, no requeue.
- Corresponding-task dependency: `aftercorr:21939`.
- Execution reports: `benchmarks/results/qfo_parameter_assessment_v1/<variant>/`.
- Short Darwin-compatible work namespaces: `qfo_benchmark/w/qpv<index>`.
- Native scoring outputs: `qfo_benchmark/scoring/parameter_v1_<index>`.
- Batch logs: `benchmarks/work/qfo_parameter_assessment_21944_<index>.log`.

All four tasks confirmed PENDING (Dependency); task 0 has raw job ID 21945.
The runner requires completed conversion, pinned converter and reference
mapping, matching count sidecar, unchanged filtered-pair content and the
fresh-native-admission record. It retains the frozen GO, EC, VGNC,
SwissTrees, TreeFam-A and FAS assessment command and input environment.
Success is explicitly pending independent score admission. Existing
execution/work/result namespaces are refused, not resumed or overwritten.

Verification: 99 focused tests passed in 0.50 seconds; batch shell syntax
passed. No parameter score has been inspected or admitted. Independent
assessment validation, paired uncertainty and the separately prespecified
CPM variants remain unfinished.

## Accounting-Stability Dependency Overrides

Before any parameter task executed, code review found that strict fresh
native-report equality and scoring-provenance checks include whole-array
`sacct` text. Capturing that text before all array members terminate could
cause a false difference when another member later completes. On September
19, the following `scontrol update` was applied separately to each index
0, 1, 2 and 3 (the commands below show index 0):

```bash
scontrol update JobId=21935_0 Dependency=afterany:21932,aftercorr:21932
scontrol update JobId=21944_0 Dependency=afterany:21939,aftercorr:21939
```

These final combined dependencies replace the original corresponding-only
dependencies in the submission sections above. Each requires the entire
upstream array to become terminal AND the corresponding task to succeed.
All eight tasks were confirmed PENDING with the combined dependencies.
An interim whole-array `afterok` setting was superseded before execution.
No frozen script/source/input changed and no running task was interrupted.
If an upstream failure leaves a dependent task permanently ineligible,
retain and audit that failure explicitly; do not bypass the success gate.

## Independent Score Admission Array 21948

Submitted 2026-09-19 at 08:49 EDT from detached executor
`benchmarks/work/publication_qfo_parameter_score_admission_v1`, commit
`d7baf35a7aa45814623debb34cdbdbb1581eff17` (pushed), passing the exact
executor path and commit to the frozen batch.

- Validator SHA-256:
  `6574b40cc97624c937eb72e2ed27c0822aff0b6ff69156815b0255055d62d81c`.
- Batch `qfo_parameter_score_admit_batch_20260919.sh` SHA-256:
  `1b551fa023ed7e1fa5094cb7393726f2962e6180108ad31493e88fa3cbb3334a`.
- Four tasks, concurrency one, 2 CPUs, 64G, 4 hours, bizon, no requeue.
- Corresponding-task dependency: `aftercorr:21944`.
- Reports: `benchmarks/work/qfo_parameter_score_admission_21948_<index>.json`.
- Logs: `benchmarks/work/qfo_parameter_score_admit_21948_<index>.log`.

Confirmed all tasks PENDING (Dependency); task 0 raw job identity is 21949.
Focused score-validation/runner/native-endpoint tests: 106 passed in 0.58
seconds; batch shell syntax passed. The validator checks complete native
output inventory, a unique successful trace, frozen command/input/helper
provenance and all six endpoint results before marking accuracy admitted.
This is not independent biological confirmation or paired uncertainty.
No parameter score is admitted merely because this array is queued.
