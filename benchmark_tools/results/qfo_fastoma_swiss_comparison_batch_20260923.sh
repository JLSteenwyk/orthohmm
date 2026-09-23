#!/bin/bash
#SBATCH --job-name=qfo_fastoma_swiss
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_fastoma_swiss_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=$ROOT/benchmarks/work/publication_corrected_swiss_comparison_v1
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == 10338e2a5046e522f2c1990e791532d3376e8982 ]]
git -C "$EXECUTOR" diff --quiet HEAD -- benchmark_tools
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
unset PYTHONPATH PYTHONHOME LD_PRELOAD LD_LIBRARY_PATH
cd "$ROOT"
exec /home/bizon/anaconda3/bin/python -B "$EXECUTOR/benchmark_tools/run_corrected_swiss_comparison.py" \
  --comparison "$ROOT/benchmark_tools/results/qfo_corrected_comparison_20260923_v6/manifest.json" \
  --comparison-sha256 6aec972b23a8ed213e5b0e2f5965bc7a16987e2d1f9a0353383ac7a8c1999c06 \
  --baseline "$ROOT/benchmark_tools/results/qfo_swiss_counts_20260917.json" \
  --protocol "$ROOT/benchmark_tools/results/QFO_SWISS_COMPARATOR_UNCERTAINTY_PROTOCOL_20260917.md" \
  --release-protocol "$ROOT/benchmark_tools/results/QFO_CORRECTED_RELEASE_PROTOCOL_20260917.md" \
  --output "$ROOT/benchmarks/work/qfo_fastoma_swiss_uncertainty_${SLURM_JOB_ID}.json"
