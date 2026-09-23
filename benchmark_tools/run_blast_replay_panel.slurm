#!/bin/bash
#SBATCH --job-name=qfo_blast_replay_diagnostic
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=180
#SBATCH --mem=900G
#SBATCH --time=02:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_blast_replay_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Pinned executor required}
COMMIT=${2:?Exact revision required}
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --exit-code HEAD -- benchmark_tools
cd "$ROOT"
unset LD_PRELOAD LD_LIBRARY_PATH LD_AUDIT PYTHONPATH
exec /home/bizon/anaconda3/bin/python "$EXECUTOR/benchmark_tools/run_blast_replay_panel.py" \
    --panel "$ROOT/benchmark_tools/results/qfo_blast_replay_panel_20260923.json"
