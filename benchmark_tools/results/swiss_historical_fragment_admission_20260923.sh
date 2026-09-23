#!/bin/bash
#SBATCH --job-name=swiss_fragment_admission
#SBATCH --nodelist=bizon
#SBATCH --dependency=afterany:22116
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/swiss_fragment_admission_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR="$ROOT/benchmarks/work/swiss_historical_fragment_admission_v1_20260923"
test "$(git -C "$EXECUTOR" rev-parse HEAD)" = "26fcd04b4887f7d667850c7c4165e7399727c332"
git -C "$EXECUTOR" diff --exit-code HEAD -- benchmark_tools
unset LD_PRELOAD LD_LIBRARY_PATH LD_AUDIT PYTHONPATH PYTHONHOME
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$EXECUTOR"
exec /home/bizon/anaconda3/bin/python -B -m benchmark_tools.verify_swiss_historical_fragments \
  --root "$ROOT" --output "$ROOT/benchmarks/work/swiss_historical_fragment_admission_20260923.json"
