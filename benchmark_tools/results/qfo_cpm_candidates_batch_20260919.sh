#!/bin/bash
#SBATCH --job-name=qfo_cpm_candidates
#SBATCH --nodelist=bizon
#SBATCH --array=0-1%1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --dependency=afterany:21962,aftercorr:21962
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_cpm_candidates_%A_%a.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Frozen executor required}
COMMIT=${2:?Exact executor commit required}
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --quiet HEAD -- benchmark_tools orthohmm
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONNOUSERSITE=1
unset PYTHONPATH PYTHONHOME LD_PRELOAD LD_LIBRARY_PATH
cd "$ROOT"
exec /usr/bin/time -v -o "$ROOT/benchmarks/work/qfo_cpm_candidates_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.time.txt" \
  /home/bizon/anaconda3/bin/python -B "$EXECUTOR/benchmark_tools/prepare_qfo_cpm_candidates.py" \
  --root "$ROOT" --index "$SLURM_ARRAY_TASK_ID"
