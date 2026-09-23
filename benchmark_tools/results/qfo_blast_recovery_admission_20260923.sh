#!/bin/bash
#SBATCH --job-name=qfo_blast_recovery_admit
#SBATCH --nodelist=bizon
#SBATCH --array=0-19%1
#SBATCH --dependency=afterany:22103,aftercorr:22103
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_blast_recovery_admit_%A_%a.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR="$ROOT/benchmarks/work/blast_recovery_admission_v1_20260923"
test "$(git -C "$EXECUTOR" rev-parse HEAD)" = "7306c5540e579854800c7a944b570ad478124485"
git -C "$EXECUTOR" diff --exit-code HEAD -- benchmark_tools
unset LD_PRELOAD LD_LIBRARY_PATH LD_AUDIT PYTHONPATH PYTHONHOME
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONHASHSEED=0
cd "$EXECUTOR"
exec /home/bizon/anaconda3/bin/python -B -m benchmark_tools.admit_blast_recovery_batch \
  --root "$ROOT" --index "$SLURM_ARRAY_TASK_ID" \
  --output "$ROOT/benchmarks/work/qfo_blast_recovery_admission_${SLURM_ARRAY_TASK_ID}_20260923.json"
