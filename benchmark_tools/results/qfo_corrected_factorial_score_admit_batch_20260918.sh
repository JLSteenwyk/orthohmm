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
[[ $# == 5 && "$COMMIT" =~ ^[0-9a-f]{40}$ && "$INDEX" =~ ^[0-7]$ ]]
[[ "$JOB" =~ ^[0-9]+$ && "$CONVERSION_JOB" =~ ^[0-9]+$ ]]
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --quiet HEAD -- benchmark_tools
CELL="p$((INDEX / 4))_c$((INDEX / 2 % 2))_r$((INDEX % 2))"
MANIFEST="$ROOT/benchmarks/results/qfo_corrected_factorial_pairs_v1/$CELL/results.json"
read -r SHA _ < <(sha256sum "$MANIFEST")
OUTPUT="$ROOT/benchmarks/work/qfo_corrected_factorial_score_admission_${INDEX}_20260918.json"
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$ROOT"
exec /home/bizon/anaconda3/bin/python "$EXECUTOR/benchmark_tools/admit_qfo_corrected_factorial_assessment.py" \
    --root "$ROOT" --index "$INDEX" --job "$JOB" --conversion-job "$CONVERSION_JOB" \
    --pairs-sha256 "$SHA" --output "$OUTPUT"
