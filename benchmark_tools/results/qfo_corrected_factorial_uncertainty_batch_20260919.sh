#!/bin/bash
#SBATCH --job-name=qfo_corrected_factorial_uncertainty
#SBATCH --partition=gpu
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_corrected_factorial_uncertainty_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Frozen analysis executor required}
COMMIT=${2:?Exact executor revision required}
[[ $# == 2 && "$COMMIT" =~ ^[0-9a-f]{40}$ ]]
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --quiet HEAD -- benchmark_tools
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$ROOT"
OUTPUT="$ROOT/benchmarks/work/qfo_corrected_factorial_uncertainty_v1"
[[ ! -e "$OUTPUT" && ! -L "$OUTPUT" ]]
ASSESSMENTS=()
for INDEX in {0..7}; do
    ADMISSION="$ROOT/benchmarks/work/qfo_corrected_factorial_score_admission_${INDEX}_20260918.json"
    [[ -f "$ADMISSION" && ! -L "$ADMISSION" ]]
    read -r SHA _ < <(sha256sum "$ADMISSION")
    [[ "$SHA" =~ ^[0-9a-f]{64}$ ]]
    ASSESSMENTS+=(--assessment "$ADMISSION" "$SHA")
done
/home/bizon/anaconda3/bin/python "$EXECUTOR/benchmark_tools/export_qfo_corrected_factorial.py" \
    "${ASSESSMENTS[@]}" --output "$OUTPUT/scores"
INVENTORY="$OUTPUT/scores/manifest.json"
read -r INVENTORY_SHA _ < <(sha256sum "$INVENTORY")
COUNTS="$OUTPUT/swiss_counts.json"
/home/bizon/anaconda3/bin/python "$EXECUTOR/benchmark_tools/audit_qfo_corrected_factorial_swiss.py" \
    --inventory "$INVENTORY" --inventory-sha256 "$INVENTORY_SHA" \
    --baseline "$ROOT/benchmark_tools/results/qfo_swiss_counts_20260917.json" --output "$COUNTS"
read -r COUNTS_SHA _ < <(sha256sum "$COUNTS")
exec /home/bizon/anaconda3/bin/python "$EXECUTOR/benchmark_tools/bootstrap_qfo_corrected_factorial.py" \
    --counts "$COUNTS" --counts-sha256 "$COUNTS_SHA" \
    --protocol "$ROOT/benchmark_tools/results/QFO_FACTORIAL_PROTOCOL_20260917.md" \
    --corrected-protocol "$ROOT/benchmark_tools/results/QFO_CORRECTED_RELEASE_PROTOCOL_20260917.md" \
    --output "$OUTPUT/swiss_bootstrap.json" --markdown "$OUTPUT/swiss_bootstrap.md"
