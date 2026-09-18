#!/bin/bash
#SBATCH --job-name=qfo_corrected_sonic
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=32
#SBATCH --mem=192G
#SBATCH --time=72:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_corrected_sonic_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Pinned executor required}
COMMIT=${2:?Exact revision required}
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --exit-code HEAD -- benchmark_tools
cd "$ROOT"
export PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
export PYTHONPYCACHEPREFIX="$ROOT/benchmarks/results/qfo_corrected_sonic_v1_pycache"
unset PYTHONPATH LD_PRELOAD LD_LIBRARY_PATH LD_AUDIT CONDA_PREFIX
exec /home/bizon/anaconda3/bin/python -B -s "$EXECUTOR/benchmark_tools/run_qfo_corrected_sonic.py" \
    --plan "$ROOT/benchmark_tools/results/qfo_corrected_sonic_commands_20260918.json"
