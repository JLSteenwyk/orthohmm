#!/bin/bash
#SBATCH --job-name=qfo_factorial_uncertainty
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_factorial_uncertainty_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Frozen analysis executor required}
COMMIT=${2:?Exact executor revision required}
INVENTORY_SHA=${3:?Reviewed eight-cell inventory SHA required}
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --quiet HEAD -- benchmark_tools
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$ROOT"
COUNTS="$ROOT/benchmarks/work/qfo_factorial_swiss_counts_20260918.json"
/home/bizon/anaconda3/bin/python "$EXECUTOR/benchmark_tools/audit_qfo_factorial_swiss.py" \
    --manifest "$ROOT/benchmark_tools/results/qfo_factorial_admission_inventory_20260918.json" \
    --manifest-sha256 "$INVENTORY_SHA" --baseline-counts "$ROOT/benchmark_tools/results/qfo_swiss_counts_20260917.json" \
    --output "$COUNTS"
HASH_RECORD=$(sha256sum "$COUNTS")
exec /home/bizon/anaconda3/bin/python "$EXECUTOR/benchmark_tools/bootstrap_qfo_factorial.py" \
    --counts "$COUNTS" --counts-sha256 "${HASH_RECORD%% *}" \
    --protocol "$ROOT/benchmark_tools/results/QFO_FACTORIAL_PROTOCOL_20260917.md" \
    --output "$ROOT/benchmarks/work/qfo_factorial_swiss_bootstrap_20260918.json" \
    --markdown "$ROOT/benchmarks/work/qfo_factorial_swiss_bootstrap_20260918.md"
