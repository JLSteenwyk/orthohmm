#!/usr/bin/env bash
#SBATCH --job-name=matched_search_graph
#SBATCH --partition=gpu
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --array=0-69%2
set -euo pipefail
readonly EXECUTOR="${1:?Pinned executor required}"
readonly COMMIT="${2:?Pinned commit required}"
readonly MANIFEST="${3:?Numeric manifest required}"
readonly INSTALL="${4:?Installation receipt required}"
readonly OUTPUT="${5:?Output root required}"
[[ $(git -C "${EXECUTOR}" rev-parse HEAD) == "${COMMIT}" ]]
git -C "${EXECUTOR}" diff --quiet HEAD -- benchmark_tools
cd "${EXECUTOR}"
exec /usr/bin/python3 -S -m benchmark_tools.run_matched_graph_cell \
    --manifest "${MANIFEST}" --installation "${INSTALL}" \
    --index "${SLURM_ARRAY_TASK_ID}" --output "${OUTPUT}/cell_${SLURM_ARRAY_TASK_ID}"
