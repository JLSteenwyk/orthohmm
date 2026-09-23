#!/bin/bash
#SBATCH --job-name=qfo_blast_recovery
#SBATCH --nodelist=bizon
#SBATCH --array=0-19%1
#SBATCH --cpus-per-task=180
#SBATCH --mem=900G
#SBATCH --time=24:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_blast_recovery_%A_%a.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR="$ROOT/benchmarks/work/blast_recovery_executor_v1_20260923"
test "$(git -C "$EXECUTOR" rev-parse HEAD)" = "1a73dc619bf50d5cc37b59f7b7820df3f4787c5d"
git -C "$EXECUTOR" diff --exit-code HEAD -- benchmark_tools
unset LD_PRELOAD LD_LIBRARY_PATH LD_AUDIT PYTHONPATH PYTHONHOME
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONHASHSEED=0
cd "$EXECUTOR"
exec /home/bizon/anaconda3/bin/python -B -m benchmark_tools.run_blast_recovery_batch \
  --root "$ROOT" --index "$SLURM_ARRAY_TASK_ID"
