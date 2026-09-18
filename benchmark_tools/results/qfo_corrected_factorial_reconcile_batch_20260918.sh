#!/bin/bash
#SBATCH --job-name=qfo_corrected_factorial_reconcile
#SBATCH --nodelist=bizon
#SBATCH --array=0-3%1
#SBATCH --cpus-per-task=32
#SBATCH --mem=192G
#SBATCH --time=48:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_corrected_factorial_reconcile_%A_%a.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Frozen executor required}
COMMIT=${2:?Exact executor revision required}
ADMISSION_JOB=${3:?Successful candidate-admission job required}
[[ $# == 3 && "$COMMIT" =~ ^[0-9a-f]{40}$ && "$ADMISSION_JOB" =~ ^[0-9]+$ ]]
[[ ${SLURM_ARRAY_TASK_ID:?Array index required} =~ ^[0-3]$ ]]
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --quiet HEAD -- benchmark_tools orthohmm
ADMISSION="$ROOT/benchmarks/work/qfo_corrected_candidate_admission_20260918.json"
read -r ADMISSION_SHA _ < <(sha256sum "$ADMISSION")
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$ROOT"
exec /home/bizon/anaconda3/bin/python "$EXECUTOR/benchmark_tools/run_qfo_corrected_factorial_cell.py" \
    --root "$ROOT" --admission "$ADMISSION" --admission-sha256 "$ADMISSION_SHA" --admission-job "$ADMISSION_JOB" \
    --environment-manifest "$ROOT/benchmark_tools/results/publication_variable_native_methods_20260916.json" \
    --index "${SLURM_ARRAY_TASK_ID:?Array index required}"
