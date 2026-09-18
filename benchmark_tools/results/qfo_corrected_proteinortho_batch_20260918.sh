#!/bin/bash
#SBATCH --job-name=qfo_corrected_proteinortho
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=32
#SBATCH --mem=192G
#SBATCH --time=72:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_corrected_proteinortho_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Pinned executor required}
COMMIT=${2:?Exact revision required}
RUNTIME_SHA=${3:?Runtime checksum required}
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --exit-code HEAD -- benchmark_tools
cd "$ROOT"
exec /home/bizon/anaconda3/bin/python "$EXECUTOR/benchmark_tools/run_qfo_corrected_proteinortho.py" run \
    --plan "$ROOT/benchmark_tools/results/qfo_corrected_proteinortho_commands_20260918.json" \
    --runtime "$ROOT/benchmark_tools/results/qfo_corrected_proteinortho_runtime_20260918.json" \
    --runtime-sha256 "$RUNTIME_SHA"
