#!/bin/bash
#SBATCH --job-name=qfo_parameter_phylogeny
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=32
#SBATCH --mem=192G
#SBATCH --time=24:00:00
#SBATCH --array=0-3%1
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_parameter_phylogeny_%A_%a.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Frozen executor required}
COMMIT=${2:?Exact executor commit required}
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --quiet HEAD -- benchmark_tools orthohmm
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
# The frozen inventory includes user-site distributions; preflight remains strict.
unset PYTHONNOUSERSITE PYTHONPATH PYTHONHOME LD_PRELOAD LD_LIBRARY_PATH
cd "$ROOT"
exec /usr/bin/time -v -o "$ROOT/benchmarks/work/qfo_parameter_phylogeny_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.time" \
    /home/bizon/anaconda3/bin/python -B "$EXECUTOR/benchmark_tools/run_qfo_parameter_phylogeny.py" \
    --root "$ROOT" --index "$SLURM_ARRAY_TASK_ID"
