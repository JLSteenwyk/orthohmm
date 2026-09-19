#!/bin/bash
#SBATCH --job-name=qfo_parameter_candidate_admit
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_parameter_candidate_admit_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Frozen admission executor required}
COMMIT=${2:?Exact admission executor commit required}
JOB=${3:?Terminal preparation job required}
MANIFEST_SHA=${4:?Completed preparation manifest SHA-256 required}
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --quiet HEAD -- benchmark_tools orthohmm
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$ROOT"
exec /home/bizon/anaconda3/bin/python "$EXECUTOR/benchmark_tools/admit_qfo_candidate_neighborhood.py" \
    --root "$ROOT" --job "$JOB" --manifest-sha256 "$MANIFEST_SHA" \
    --output "$ROOT/benchmarks/work/qfo_parameter_candidate_admission_${SLURM_JOB_ID}.json"
