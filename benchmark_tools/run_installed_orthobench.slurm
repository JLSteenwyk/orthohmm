#!/usr/bin/env bash
#SBATCH --job-name=installed_ob_repro
#SBATCH --partition=gpu
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=24:00:00
set -euo pipefail
readonly EXECUTOR="${1:?Pinned executor worktree required}"
readonly COMMIT="${2:?Pinned commit required}"
readonly PLAN="${3:?Prepared plan required}"
readonly PLAN_SHA="${4:?Prepared plan SHA-256 required}"
[[ $(sha256sum "${PLAN}") == "${PLAN_SHA}  ${PLAN}" ]]
[[ $(git -C "${EXECUTOR}" rev-parse HEAD) == "${COMMIT}" ]]
git -C "${EXECUTOR}" diff --quiet HEAD -- benchmark_tools
cd "${EXECUTOR}"
exec /usr/bin/python3 -S -m benchmark_tools.run_installed_orthobench --run "${PLAN}"
