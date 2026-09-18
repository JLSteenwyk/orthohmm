#!/bin/bash
#SBATCH --job-name=qfo_corrected_group_pairs
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_corrected_group_pairs_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Frozen executor required}
COMMIT=${2:?Exact revision required}
INDEX=${3:?R-off factorial index required}
ADMISSION=${4:?Corrected candidate admission required}
SHA=${5:?Reviewed admission SHA-256 required}
JOB=${6:?Candidate admission job required}
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --quiet HEAD -- benchmark_tools orthohmm qfo_benchmark/og_to_pairwise.py
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$ROOT"
exec /home/bizon/anaconda3/bin/python "$EXECUTOR/benchmark_tools/prepare_qfo_corrected_group_pairs.py" \
    --root "$ROOT" --index "$INDEX" --admission "$ADMISSION" \
    --admission-sha256 "$SHA" --admission-job "$JOB"
