# Corrected Replay Preparation Queued

Job **21722** is scheduler-confirmed pending `afterok:21720` on bizon,
with 2 CPUs, 64 GiB and a four-hour limit, without requeue. It prepares a
command manifest only. It does not launch replay, reconciliation or scoring.

Dependency chain:

1. Corrected native OrthoHMM inference: array task **21706_0**, raw job 21707.
2. Independent native evidence admission: **21720**, after inference ends.
3. Corrected replay command preparation: **21722**, only if admission succeeds.

The detached preparation/driver executor is
`benchmarks/work/publication_qfo_corrected_replay_v1`, revision
`188fde21860a70da55fee1358177485c970a11a3`.
Batch: `qfo_corrected_replay_prepare_batch_20260918.sh`.
Input: `benchmarks/work/qfo_corrected_high_sensitivity_admission_20260918.json`.
Expected plan: `benchmarks/work/qfo_corrected_replay_commands_20260918.json`.
Reserved fresh output: `benchmarks/results/qfo_corrected_checked_replay_v1`.
Log: `benchmarks/work/qfo_corrected_replay_prepare_21722.log`.

The batch captures the actual completed admission file SHA-256; the builder
then checks the admitted corrected FASTAs/checkpoint, fixed primary plan and
frozen launcher/runtime. Existing destinations are refused. A failed native
admission leaves this job in dependency failure rather than producing a plan
from partial inputs. Thirty-nine focused preparation/worker/driver tests and
batch shell syntax checks passed before submission.

## Remaining Execution Gates

Neither an admission report nor corrected replay plan is available at
submission. Inspect successful preparation and the actual plan hash, preserve
the manifest, then freeze the 32-CPU checked replay launch. Do not run the
manifest's underlying replay command directly: it must pass through
`run_qfo_corrected_replay.py` and the four checked native graph boundaries.

Independent replay admission must verify retained per-stage partitions, source,
runtime, command, full corrected gene coverage and native partition comparison
before factorial candidate generation. Original-release replay artifacts cannot
supply these corrected stages. Nonequivalence must be preserved rather than
resolved through retries or score transfer. Resource costs remain incremental
shared-host evidence, separate from DGX end-to-end timings.
