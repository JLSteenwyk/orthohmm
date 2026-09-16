#!/usr/bin/env bash
#SBATCH --job-name=orthobench_factorial
#SBATCH --partition=gpu
#SBATCH --nodelist=bizon
#SBATCH --array=0-3%2
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=24:00:00
set -euo pipefail

: "${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID required}"
readonly ROOT="${1:?Repository root required}"
readonly EXECUTOR="${2:?Pinned executor required}"
readonly COMMIT="${3:?Pinned executor commit required}"
shift 3
[[ $# == 0 || ( $# == 1 && $1 == --check-only ) ]] || exit 2
[[ $(git -C "${EXECUTOR}" rev-parse HEAD) == "${COMMIT}" ]]
git -C "${EXECUTOR}" diff --quiet HEAD -- benchmark_tools
exec /home/bizon/anaconda3/bin/python "${EXECUTOR}/benchmark_tools/run_orthobench_factorial_cell.py" \
    --manifest "${ROOT}/benchmark_tools/results/orthobench_factorial_prepared_20260916.json" \
    --manifest-sha256 5c325f4d77865e0c7571fe4bb4df0be0959977a49e1d22169f192518700f9382 \
    --environment-manifest "${ROOT}/benchmark_tools/results/publication_variable_native_methods_20260916.json" \
    --environment-sha256 bf677728cf9312edd9755d9ac1ceefbce3abb8c3848d43413631dc00ca00eb2f \
    --preparation-job 21161 --index "${SLURM_ARRAY_TASK_ID}" "$@"
