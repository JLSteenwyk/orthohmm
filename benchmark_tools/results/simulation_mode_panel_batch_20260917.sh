#!/bin/bash
#SBATCH --job-name=simulation_mode_panel
#SBATCH --array=0-69%2
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/simulation_mode_panel_%A_%a.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Frozen executor required}
COMMIT=${2:?Frozen revision required}
MANIFEST=${3:?Frozen panel manifest required}
HASH=${4:?Frozen panel SHA256 required}
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --quiet HEAD -- benchmark_tools
cd "$ROOT"
exec /home/bizon/anaconda3/bin/python "$EXECUTOR/benchmark_tools/run_simulation_mode_panel.py" \
    --root "$ROOT" --manifest "$MANIFEST" --manifest-sha256 "$HASH" --index "$SLURM_ARRAY_TASK_ID"
