#!/bin/bash
#SBATCH --job-name=qfo_corrected_legacy_blast
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=180
#SBATCH --mem=900G
#SBATCH --time=14-00:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_corrected_blast_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Pinned executor required}
COMMIT=${2:?Exact revision required}
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --exit-code HEAD -- benchmark_tools
cd "$ROOT"
unset LD_PRELOAD LD_LIBRARY_PATH LD_AUDIT PYTHONPATH
exec /home/bizon/anaconda3/bin/python "$EXECUTOR/benchmark_tools/run_qfo_corrected_blast.py" \
    --plan "$ROOT/benchmark_tools/results/qfo_corrected_orthomcl_prepared_20260918.json" \
    --runtime "$ROOT/benchmark_tools/results/qfo_corrected_legacy_blast_runtime_20260918.json"
