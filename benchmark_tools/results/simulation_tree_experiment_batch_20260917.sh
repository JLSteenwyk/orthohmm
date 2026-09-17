#!/bin/bash
#SBATCH --job-name=simulation_tree_experiment
#SBATCH --array=0-209%2
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/simulation_tree_experiment_%A_%a.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Frozen executor required}
COMMIT=${2:?Frozen executor revision required}
ADMISSION_SHA=${3:?Reviewed completed mode-panel admission SHA256 required}
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --quiet HEAD -- benchmark_tools
[[ ${SLURM_ARRAY_TASK_ID:?} =~ ^[0-9]+$ ]]
[[ $ADMISSION_SHA =~ ^[a-f0-9]{64}$ ]]
cd "$ROOT"
exec /home/bizon/anaconda3/bin/python "$EXECUTOR/benchmark_tools/run_simulation_tree_experiment.py" \
    --root "$ROOT" --index "$SLURM_ARRAY_TASK_ID" --mode-admission-sha256 "$ADMISSION_SHA" \
    --output "$ROOT/benchmarks/results/simulation_tree_experiments_v1/cell_${SLURM_ARRAY_TASK_ID}"
