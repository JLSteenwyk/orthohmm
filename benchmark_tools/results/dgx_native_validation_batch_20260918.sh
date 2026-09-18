#!/bin/bash
#SBATCH --job-name=dgx_native_output_audit
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/dgx_native_validation_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Frozen executor required}
COMMIT=${2:?Exact revision required}
[[ $# == 2 && "$COMMIT" =~ ^[0-9a-f]{40}$ ]]
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --exit-code HEAD -- benchmark_tools orthohmm
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$ROOT"
exec /home/bizon/anaconda3/bin/python "$EXECUTOR/benchmark_tools/validate_dgx_native_panel.py" \
    --archive "$ROOT/benchmarks/work/dgx_native_archive_20260918" \
    --results "$ROOT/benchmark_tools/results" \
    --inventory "$ROOT/benchmarks/work/dgx_native_archive_inventory_20260918.json" \
    --output "$ROOT/benchmarks/work/dgx_native_output_validation_20260918.json"
