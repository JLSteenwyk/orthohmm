#!/bin/bash
#SBATCH --job-name=qfo_corrected_replay_prepare
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_corrected_replay_prepare_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Pinned preparation executor required}
COMMIT=${2:?Exact executor revision required}
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --exit-code HEAD -- benchmark_tools orthohmm
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$ROOT"
ADMISSION="$ROOT/benchmarks/work/qfo_corrected_high_sensitivity_admission_20260918.json"
# Submit only after successful independent admission. The builder rechecks
# admitted inputs and freezes this content hash; this batch does not run replay.
HASH_RECORD=$(sha256sum "$ADMISSION")
ADMISSION_SHA=${HASH_RECORD%% *}
exec /home/bizon/anaconda3/bin/python "$EXECUTOR/benchmark_tools/prepare_qfo_corrected_replay.py" \
    --root "$ROOT" --admission "$ADMISSION" --admission-sha256 "$ADMISSION_SHA" \
    --output-root "$ROOT/benchmarks/results/qfo_corrected_checked_replay_v1" \
    --manifest "$ROOT/benchmarks/work/qfo_corrected_replay_commands_20260918.json"
