#!/bin/bash
#SBATCH --job-name=qfo_factorial_reconcile
#SBATCH --nodelist=bizon
#SBATCH --array=0-3%1
#SBATCH --cpus-per-task=32
#SBATCH --mem=192G
#SBATCH --time=48:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_factorial_reconcile_%A_%a.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Frozen executor required}
COMMIT=${2:?Exact executor revision required}
MANIFEST=${3:?Prepared manifest required}
MANIFEST_SHA=${4:?Prepared manifest SHA256 required}
PREPARATION_JOB=${5:?Completed preparation job required}
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --quiet HEAD -- benchmark_tools
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$ROOT"
exec /home/bizon/anaconda3/bin/python "$EXECUTOR/benchmark_tools/run_qfo_factorial_cell.py" \
    --manifest "$MANIFEST" --manifest-sha256 "$MANIFEST_SHA" \
    --environment-manifest "$ROOT/benchmark_tools/results/publication_variable_native_methods_20260916.json" \
    --preparation-job "$PREPARATION_JOB" --index "${SLURM_ARRAY_TASK_ID:?Array index required}"
