#!/bin/bash
#SBATCH --job-name=qfo_corrected_score_admit
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_corrected_score_admit_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Pinned executor required}
COMMIT=${2:?Exact revision required}
METHOD=${3:?Method required}
JOB=${4:?Assessment job required}
CONVERSION_JOB=${5:?Conversion job required}
SHA=${6:?Reviewed pair-manifest checksum required}
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --exit-code HEAD -- benchmark_tools
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$ROOT"
exec /home/bizon/anaconda3/bin/python "$EXECUTOR/benchmark_tools/admit_qfo_corrected_comparator_assessment.py" \
    --root "$ROOT" --method "$METHOD" --job "$JOB" --conversion-job "$CONVERSION_JOB" \
    --pairs-sha256 "$SHA" --output "$ROOT/benchmarks/work/qfo_corrected_${METHOD}_assessment_admission_20260918.json"
