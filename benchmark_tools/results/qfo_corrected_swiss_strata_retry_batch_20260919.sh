#!/bin/bash
#SBATCH --job-name=qfo_corrected_swiss_strata_v2
#SBATCH --partition=gpu
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --no-requeue
#SBATCH --dependency=afterok:21894:21736_0
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_corrected_swiss_strata_retry_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Frozen analysis executor required}
COMMIT=${2:?Exact executor revision required}
[[ $# == 2 && "$COMMIT" =~ ^[0-9a-f]{40}$ ]]
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --quiet HEAD -- benchmark_tools
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$ROOT"
INVENTORY="$ROOT/benchmarks/work/qfo_corrected_factorial_uncertainty_v1/scores/manifest.json"
ORTHOFINDER="$ROOT/benchmarks/work/qfo_corrected_orthofinder_full_assessment_admission_20260918.json"
OUTPUT="$ROOT/benchmarks/work/qfo_corrected_swiss_primary_strata_v2.json"
[[ ! -e "$OUTPUT" && ! -L "$OUTPUT" ]]
[[ -f "$INVENTORY" && ! -L "$INVENTORY" && -f "$ORTHOFINDER" && ! -L "$ORTHOFINDER" ]]
read -r INVENTORY_SHA _ < <(sha256sum "$INVENTORY")
read -r ORTHOFINDER_SHA _ < <(sha256sum "$ORTHOFINDER")
[[ "$INVENTORY_SHA" =~ ^[0-9a-f]{64}$ && "$ORTHOFINDER_SHA" =~ ^[0-9a-f]{64}$ ]]
/home/bizon/anaconda3/bin/python -c 'import json,platform,numpy; assert platform.python_version() == "3.10.13"; assert numpy.__version__ == "2.2.6"; print(json.dumps({"python":platform.python_version(),"numpy":numpy.__version__}),flush=True)'
exec /home/bizon/anaconda3/bin/python "$EXECUTOR/benchmark_tools/run_corrected_swiss_strata.py" \
    --inventory "$INVENTORY" --inventory-sha256 "$INVENTORY_SHA" \
    --orthofinder-admission "$ORTHOFINDER" --orthofinder-admission-sha256 "$ORTHOFINDER_SHA" \
    --baseline "$ROOT/benchmark_tools/results/qfo_swiss_counts_20260917.json" \
    --strata "$ROOT/benchmark_tools/results/corrected_swiss_sequence_strata_20260918.json" \
    --protocol "$ROOT/benchmark_tools/results/CORRECTED_SWISS_SEQUENCE_STRATA_PROTOCOL_20260918.md" \
    --output "$OUTPUT"
