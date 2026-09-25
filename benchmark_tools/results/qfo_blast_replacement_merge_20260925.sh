#!/bin/bash
#SBATCH --job-name=qfo_blast_replacement_merge
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --dependency=afterok:22161:22105_19
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_blast_replacement_merge_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR="$ROOT/benchmarks/work/blast_replacement_merge_v1_20260925"
test "$(git -C "$EXECUTOR" rev-parse HEAD)" = f6ab36db0cbed151ac2b46581345cd187d514e79
test -z "$(git -C "$EXECUTOR" status --porcelain)"
[[ "$SLURM_JOB_NODELIST" == bizon && "$SLURM_CPUS_PER_TASK" == 2 ]]
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
unset PYTHONPATH PYTHONHOME LD_PRELOAD LD_LIBRARY_PATH LD_AUDIT
cd "$EXECUTOR"
exec /home/bizon/anaconda3/bin/python -B -m benchmark_tools.run_blast_recovery_merge \
  --root "$ROOT" --replacement --output "$ROOT/benchmarks/results/qfo_blast_replacement_merge_v1"
