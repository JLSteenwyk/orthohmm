#!/bin/bash
#SBATCH --job-name=qfo_blast_replacement_search_admit
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --dependency=afterok:22162
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_blast_replacement_search_admit_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Frozen executor required}
COMMIT=${2:?Exact executor commit required}
test "$(git -C "$EXECUTOR" rev-parse HEAD)" = "$COMMIT"
test -z "$(git -C "$EXECUTOR" status --porcelain)"
[[ "$SLURM_JOB_NODELIST" == bizon && "$SLURM_CPUS_PER_TASK" == 2 ]]
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
unset PYTHONPATH PYTHONHOME LD_PRELOAD LD_LIBRARY_PATH LD_AUDIT
cd "$EXECUTOR"
exec /home/bizon/anaconda3/bin/python -B -m benchmark_tools.admit_blast_recovery_search \
  --root "$ROOT" --replacement --output "$ROOT/benchmarks/results/qfo_blast_replacement_search_admission_v1"
