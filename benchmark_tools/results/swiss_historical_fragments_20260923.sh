#!/bin/bash
#SBATCH --job-name=swiss_historical_fragments
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/swiss_historical_fragments_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR="$ROOT/benchmarks/work/swiss_historical_fragments_v1_20260923"
test "$(git -C "$EXECUTOR" rev-parse HEAD)" = "df8a9c09a8c7a1f4476d1d6a7a9c6ccb538e69b1"
git -C "$EXECUTOR" diff --exit-code HEAD -- benchmark_tools
unset LD_PRELOAD LD_LIBRARY_PATH LD_AUDIT PYTHONPATH PYTHONHOME
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$EXECUTOR"
exec /home/bizon/anaconda3/bin/python -B -m benchmark_tools.collect_swiss_historical_fragments \
  --root "$ROOT" --output "$ROOT/benchmarks/work/swiss_historical_fragment_panel_20260923"
