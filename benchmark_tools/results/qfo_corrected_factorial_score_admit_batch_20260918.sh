#!/bin/bash
#SBATCH --job-name=qfo_corrected_factorial_score_admit
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_corrected_factorial_score_admit_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Frozen admission executor required}
COMMIT=${2:?Exact revision required}
INDEX=${3:?Corrected cell index required}
JOB=${4:?Completed assessment job required}
CONVERSION_JOB=${5:?Completed conversion job required}
SHA=${6:?Reviewed pair manifest SHA-256 required}
OUTPUT=${7:?Fresh admission output required}
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --quiet HEAD -- benchmark_tools
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$ROOT"
exec /home/bizon/anaconda3/bin/python "$EXECUTOR/benchmark_tools/admit_qfo_corrected_factorial_assessment.py" \
    --root "$ROOT" --index "$INDEX" --job "$JOB" --conversion-job "$CONVERSION_JOB" \
    --pairs-sha256 "$SHA" --output "$OUTPUT"
