#!/bin/bash
#SBATCH --job-name=ob_search_decisions
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/ob_search_decisions_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Frozen diagnostic executor required}
COMMIT=${2:?Exact executor revision required}
[[ $# == 2 && "$COMMIT" =~ ^[0-9a-f]{40}$ ]]
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --exit-code HEAD -- benchmark_tools
export PYTHONHASHSEED=0 OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=4
exec /home/bizon/anaconda3/bin/python "$EXECUTOR/benchmark_tools/trace_ob_search_decisions.py" \
    --root "$ROOT" --core "$ROOT/benchmarks/work/publication_method_native_v2" \
    --output "$ROOT/benchmarks/results/ob_search_decisions_v1" --threads 4
