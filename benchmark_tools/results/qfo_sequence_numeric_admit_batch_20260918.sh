#!/bin/bash
#SBATCH --job-name=qfo_corrected_numeric_admit
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=2
#SBATCH --mem=192G
#SBATCH --time=7-00:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_sequence_numeric_admit_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Frozen validator required}
COMMIT=${2:?Exact revision required}
JOB=${3:?Conversion job required}
[[ $# == 3 && "$COMMIT" =~ ^[0-9a-f]{40}$ && "$JOB" =~ ^[0-9]+$ ]]
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --exit-code HEAD -- benchmark_tools orthohmm
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$ROOT"
exec /home/bizon/anaconda3/bin/python "$EXECUTOR/benchmark_tools/admit_qfo_sequence_numeric.py" \
    --root "$ROOT" --job "$JOB" \
    --output "$ROOT/benchmarks/work/qfo_sequence_numeric_admission_20260918"
