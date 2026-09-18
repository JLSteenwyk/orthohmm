#!/bin/bash
#SBATCH --job-name=qfo_corrected_hit_coverage
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=2
#SBATCH --mem=192G
#SBATCH --time=7-00:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_search_coverage_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Frozen executor required}
COMMIT=${2:?Exact revision required}
HMM_JOB=${3:?HMM admission job required}
NUMERIC_JOB=${4:?Numeric admission job required}
[[ $# == 4 && "$COMMIT" =~ ^[0-9a-f]{40}$ && "$HMM_JOB" =~ ^[0-9]+$ && "$NUMERIC_JOB" =~ ^[0-9]+$ ]]
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --exit-code HEAD -- benchmark_tools orthohmm
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$ROOT"
exec /home/bizon/anaconda3/bin/python "$EXECUTOR/benchmark_tools/compare_qfo_search_coverage.py" \
    --root "$ROOT" --hmm-job "$HMM_JOB" --numeric-job "$NUMERIC_JOB" \
    --output "$ROOT/benchmarks/work/qfo_corrected_hit_coverage_20260918"
