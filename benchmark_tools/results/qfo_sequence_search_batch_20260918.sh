#!/bin/bash
#SBATCH --job-name=qfo_corrected_sequence_search
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=32
#SBATCH --mem=192G
#SBATCH --time=7-00:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_sequence_search_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Frozen executor required}
COMMIT=${2:?Exact revision required}
[[ $# == 2 && "$COMMIT" =~ ^[0-9a-f]{40}$ ]]
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --exit-code HEAD -- benchmark_tools orthohmm
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$ROOT"
exec /home/bizon/anaconda3/bin/python "$EXECUTOR/benchmark_tools/run_qfo_sequence_search_control.py" \
    --manifest "$ROOT/benchmarks/work/qfo_sequence_search_control_v1/manifest.json"
