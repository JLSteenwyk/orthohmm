#!/bin/bash
#SBATCH --job-name=qfo_label_recovery
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_label_recovery_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Frozen auditor required}
COMMIT=${2:?Frozen auditor revision required}
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --quiet HEAD -- benchmark_tools
cd "$ROOT"
exec /home/bizon/anaconda3/bin/python "$EXECUTOR/benchmark_tools/admit_qfo_checked_full_replay.py" \
    --root "$ROOT" --output "$ROOT/benchmarks/results/qfo_checked_full_replay_v2_recovered_admission_v1" \
    --report-sha256 6a6b682e2e86566e59ddf8748f6eeb475eab976ee35c988df73bf1b75d4c100a \
    --run-version v2 --recover-label-failure
