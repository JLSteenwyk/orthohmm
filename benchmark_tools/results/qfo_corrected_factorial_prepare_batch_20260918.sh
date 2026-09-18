#!/bin/bash
#SBATCH --job-name=qfo_corrected_factorial_prepare
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_corrected_factorial_prepare_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Frozen preparation executor required}
COMMIT=${2:?Exact executor revision required}
ADMISSION=${3:?Corrected replay admission required}
ADMISSION_SHA=${4:?Reviewed admission SHA-256 required}
ADMISSION_JOB=${5:?Successful independent admission job required}
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --quiet HEAD -- benchmark_tools orthohmm
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$ROOT"
exec /home/bizon/anaconda3/bin/python "$EXECUTOR/benchmark_tools/prepare_qfo_corrected_factorial.py" \
    --root "$ROOT" --admission "$ADMISSION" --admission-sha256 "$ADMISSION_SHA" \
    --admission-job "$ADMISSION_JOB" --output "$ROOT/benchmarks/results/qfo_corrected_factorial_v1"
