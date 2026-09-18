#!/bin/bash
#SBATCH --job-name=qfo_corrected_factorial_admit
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_corrected_factorial_admit_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Frozen admission executor required}
COMMIT=${2:?Exact revision required}
INDEX=${3:?Reconciliation index required}
JOB=${4:?Raw reconciliation job ID required}
ADMISSION=${5:?Corrected candidate admission required}
SHA=${6:?Candidate admission SHA-256 required}
ADMISSION_JOB=${7:?Candidate admission job required}
OUTPUT=${8:?Fresh result path required}
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --quiet HEAD -- benchmark_tools orthohmm
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$ROOT"
exec /home/bizon/anaconda3/bin/python "$EXECUTOR/benchmark_tools/admit_qfo_corrected_factorial_cell.py" \
    --root "$ROOT" --index "$INDEX" --job "$JOB" --admission "$ADMISSION" \
    --admission-sha256 "$SHA" --admission-job "$ADMISSION_JOB" --output "$OUTPUT"
