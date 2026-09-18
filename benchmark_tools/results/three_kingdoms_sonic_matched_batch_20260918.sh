#!/bin/bash
#SBATCH --job-name=tk_sonic_matched
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=32
#SBATCH --mem=192G
#SBATCH --time=72:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/three_kingdoms_sonic_matched_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Pinned executor required}
COMMIT=${2:?Exact revision required}
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --exit-code HEAD -- benchmark_tools
cd "$ROOT"
export PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
export PYTHONPYCACHEPREFIX="$ROOT/benchmarks/results/three_kingdoms_sonic_matched_v1_pycache"
unset PYTHONPATH PYTHONHOME LD_PRELOAD LD_LIBRARY_PATH LD_AUDIT CONDA_PREFIX
exec /home/bizon/anaconda3/bin/python -B -s "$EXECUTOR/benchmark_tools/run_three_kingdoms_matched_sonic.py" run \
    --plan "$ROOT/benchmark_tools/results/three_kingdoms_sonic_matched_commands_20260918.json" \
    --sha256 5c3ce4608b54e8e550f0ea855754879211b842f0175b2e8e864cf4b5c2998283
