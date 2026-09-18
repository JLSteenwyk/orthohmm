#!/bin/bash
#SBATCH --job-name=qfo_sequence_pairs
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_sequence_pairs_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Frozen executor required}
COMMIT=${2:?Exact revision required}
ADMISSION=${3:?Admission report required}
SHA=${4:?Admission hash required}
JOB=${5:?Admission job required}
VARIANT=${6:?Variant required}
[[ $# == 6 && "$COMMIT" =~ ^[0-9a-f]{40}$ && "$SHA" =~ ^[0-9a-f]{64}$ && "$JOB" =~ ^[0-9]+$ ]]
[[ "$VARIANT" == all_hits || "$VARIANT" == top100 ]]
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --exit-code HEAD -- benchmark_tools qfo_benchmark orthohmm
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$ROOT"
exec /home/bizon/anaconda3/bin/python "$EXECUTOR/benchmark_tools/prepare_qfo_sequence_pairs.py" \
    --root "$ROOT" --admission "$ADMISSION" --admission-sha256 "$SHA" \
    --admission-job "$JOB" --variant "$VARIANT"
