#!/bin/bash
#SBATCH --job-name=ob_initial_edge_trace
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/ob_initial_edge_trace_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Executor required}
COMMIT=${2:?Revision required}
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --exit-code HEAD -- benchmark_tools
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONHASHSEED=0
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
unset PYTHONPATH PYTHONHOME LD_PRELOAD LD_LIBRARY_PATH LD_AUDIT
exec /home/bizon/anaconda3/bin/python -B -s "$EXECUTOR/benchmark_tools/trace_ob_initial_edges.py" \
  --repo "$ROOT" --core "$ROOT/benchmarks/work/publication_method_native_v2" \
  --output "$ROOT/benchmarks/work/ob_initial_edge_trace_20260918"
