#!/bin/bash
#SBATCH --job-name=qfo_cpm_checkpoint_admit
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --no-requeue
#SBATCH --dependency=afterok:22154
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_cpm_checkpoint_admit_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Frozen admission executor required}
COMMIT=${2:?Exact admission executor commit required}
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --quiet HEAD -- benchmark_tools orthohmm
export PYTHONHASHSEED=0 PYTHONFAULTHANDLER=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONNOUSERSITE=1
unset PYTHONPATH PYTHONHOME LD_PRELOAD LD_LIBRARY_PATH
cd "$ROOT"
exec /usr/bin/time -v -o "$ROOT/benchmarks/work/qfo_cpm_checkpoint_admit_${SLURM_JOB_ID}.time.txt" \
  /home/bizon/anaconda3/bin/python -B "$EXECUTOR/benchmark_tools/admit_cpm_checkpoint_recovery.py" \
  --root "$ROOT" --job 22154 --commit 0637916c14a81e7a3b31f5aeb5fadfc52b103079 \
  --output "$ROOT/benchmarks/results/qfo_cpm_checkpoint_recovery_admission_v1"
