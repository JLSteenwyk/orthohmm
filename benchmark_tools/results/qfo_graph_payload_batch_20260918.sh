#!/bin/bash
#SBATCH --job-name=qfo_graph_payload
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_graph_payload_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Frozen executor required}
COMMIT=${2:?Exact revision required}
ADMISSION_JOB=${3:?Numeric admission job required}
[[ $# == 3 && "$COMMIT" =~ ^[0-9a-f]{40}$ && "$ADMISSION_JOB" =~ ^[0-9]+$ ]]
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --exit-code HEAD -- benchmark_tools orthohmm
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$ROOT"
exec /home/bizon/anaconda3/bin/python "$EXECUTOR/benchmark_tools/plan_qfo_sequence_graph_memory.py" \
    --root "$ROOT" --admission-job "$ADMISSION_JOB" \
    --output "$ROOT/benchmarks/work/qfo_graph_payload_20260918.json"
