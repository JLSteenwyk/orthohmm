#!/bin/bash
#SBATCH --job-name=qfo_sequence_graph
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=32
#SBATCH --time=7-00:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_sequence_graph_%j.log
# Supply --mem explicitly at submission to match the reviewed per-variant plan.
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Frozen executor required}
COMMIT=${2:?Exact revision required}
PLAN=${3:?Plan required}
SHA=${4:?Plan hash required}
VARIANT=${5:?Variant required}
[[ $# == 5 && "$COMMIT" =~ ^[0-9a-f]{40}$ && "$SHA" =~ ^[0-9a-f]{64}$ ]]
[[ "$VARIANT" == all_hits || "$VARIANT" == top100 ]]
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --exit-code HEAD -- benchmark_tools orthohmm
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$ROOT"
exec /home/bizon/anaconda3/bin/python "$EXECUTOR/benchmark_tools/run_qfo_sequence_graph.py" \
    --root "$ROOT" --plan "$PLAN" --plan-sha256 "$SHA" --variant "$VARIANT"
