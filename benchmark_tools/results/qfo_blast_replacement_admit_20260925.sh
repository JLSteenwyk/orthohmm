#!/bin/bash
#SBATCH --job-name=qfo_blast_replacement_admit
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --no-requeue
#SBATCH --dependency=afterok:22160_14
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_blast_replacement_admit_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?frozen executor required}
COMMIT=${2:?frozen commit required}
test "$(git -C "$EXECUTOR" rev-parse HEAD)" = "$COMMIT"
test -z "$(git -C "$EXECUTOR" status --porcelain)"
unset LD_PRELOAD LD_LIBRARY_PATH LD_AUDIT PYTHONPATH PYTHONHOME
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONHASHSEED=0
cd "$EXECUTOR"
exec /home/bizon/anaconda3/bin/python -B -m benchmark_tools.admit_blast_recovery_batch \
  --root "$ROOT" --index 14 --replacement \
  --output "$ROOT/benchmarks/work/qfo_blast_recovery_admission_14_20260923.json"
