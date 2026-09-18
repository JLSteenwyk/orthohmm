#!/bin/bash
#SBATCH --job-name=tk_sonic_assess
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/three_kingdoms_sonic_assess_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Pinned executor required}
COMMIT=${2:?Exact revision required}
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --exit-code HEAD -- benchmark_tools
export PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
unset PYTHONPATH PYTHONHOME LD_PRELOAD LD_LIBRARY_PATH LD_AUDIT CONDA_PREFIX
exec /home/bizon/anaconda3/bin/python -B -s "$EXECUTOR/benchmark_tools/assess_three_kingdoms_matched_sonic.py" \
  --repo "$ROOT" --output "$ROOT/benchmarks/work/three_kingdoms_sonic_matched_assessment_20260918"
