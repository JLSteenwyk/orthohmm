#!/bin/bash
#SBATCH --job-name=qfo_corrected_checked_replay
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=32
#SBATCH --mem=192G
#SBATCH --time=24:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_corrected_checked_replay_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Frozen replay executor required}
COMMIT=${2:?Exact executor revision required}
PLAN=${3:?Admitted corrected replay plan required}
PLAN_SHA=${4:?Reviewed plan SHA-256 required}
[[ $# == 4 && "$COMMIT" =~ ^[0-9a-f]{40}$ && "$PLAN_SHA" =~ ^[0-9a-f]{64}$ ]]
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --quiet HEAD -- benchmark_tools orthohmm
HASH_RECORD=$(sha256sum "$PLAN")
[[ ${HASH_RECORD%% *} == "$PLAN_SHA" ]]
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$ROOT"
exec /home/bizon/anaconda3/bin/python "$EXECUTOR/benchmark_tools/run_qfo_corrected_replay.py" \
    --root "$ROOT" --plan "$PLAN" --plan-sha256 "$PLAN_SHA"
