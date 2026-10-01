#!/bin/bash
#SBATCH --job-name=qfo_private_cpm_parameter_uncertainty
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_private_cpm_parameter_uncertainty_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Frozen parameter-integration executor required}
COMMIT=${2:?Exact integration revision required}
PROTOCOL_SHA=${3:?Reviewed integration protocol required}
INVENTORY=${4:?Fresh reviewed seven-arm inventory required}
INVENTORY_SHA=${5:?Reviewed inventory SHA256 required}
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --quiet HEAD -- benchmark_tools orthohmm
export PYTHONHASHSEED=0 PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
unset PYTHONPATH PYTHONHOME LD_PRELOAD LD_LIBRARY_PATH LD_AUDIT
cd "$ROOT"
exec /home/bizon/anaconda3/bin/python -B "$EXECUTOR/benchmark_tools/run_private_cpm_parameter_uncertainty.py" \
  --root "$ROOT" --protocol-sha256 "$PROTOCOL_SHA" \
  --inventory "$INVENTORY" --inventory-sha256 "$INVENTORY_SHA" \
  --baseline "$ROOT/benchmark_tools/results/qfo_swiss_counts_20260917.json" \
  --output "$ROOT/benchmarks/work/qfo_private_cpm_parameter_uncertainty_20261001.json"
