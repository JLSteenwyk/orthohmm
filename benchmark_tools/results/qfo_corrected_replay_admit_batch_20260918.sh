#!/bin/bash
#SBATCH --job-name=qfo_corrected_replay_admit
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_corrected_replay_admit_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Frozen admission executor required}
COMMIT=${2:?Exact executor revision required}
JOB=${3:?Replay job required}
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --exit-code HEAD -- benchmark_tools orthohmm
PLAN="$ROOT/benchmarks/work/qfo_corrected_replay_commands_20260918.json"
REPORT="$ROOT/benchmarks/results/qfo_corrected_checked_replay_v1/results.json"
read -r PLAN_SHA _ < <(sha256sum "$PLAN")
read -r REPORT_SHA _ < <(sha256sum "$REPORT")
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$ROOT"
exec /home/bizon/anaconda3/bin/python "$EXECUTOR/benchmark_tools/admit_qfo_corrected_replay.py" \
    --root "$ROOT" --plan "$PLAN" --plan-sha256 "$PLAN_SHA" --report-sha256 "$REPORT_SHA" \
    --job "$JOB" --output "$ROOT/benchmarks/work/qfo_corrected_replay_admission_20260918.json"
