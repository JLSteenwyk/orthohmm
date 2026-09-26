#!/usr/bin/env bash
#SBATCH --job-name=search_sensitivity
#SBATCH --partition=gpu
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --array=0-69%2
set -euo pipefail
readonly EXECUTOR="${1:?Pinned executor worktree required}"
readonly COMMIT="${2:?Pinned commit required}"
readonly PLAN="${3:?Frozen plan required}"
readonly OUTPUT="${4:?Output root required}"
[[ $(git -C "${EXECUTOR}" rev-parse HEAD) == "${COMMIT}" ]]
git -C "${EXECUTOR}" diff --quiet HEAD -- benchmark_tools
cd "${EXECUTOR}"
exec /usr/bin/python3 -S -m benchmark_tools.run_search_sensitivity_cell \
    --plan "${PLAN}" --index "${SLURM_ARRAY_TASK_ID}" \
    --output "${OUTPUT}/cell_${SLURM_ARRAY_TASK_ID}"
