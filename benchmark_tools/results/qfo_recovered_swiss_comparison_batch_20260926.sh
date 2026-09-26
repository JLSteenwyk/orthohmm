#!/bin/bash
#SBATCH --job-name=qfo_recovered_swiss
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_recovered_swiss_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=$ROOT/benchmarks/work/publication_recovered_swiss_v1_20260923
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == 0b8efa9c6b20b5bdec41f798c136b54132e7db37 ]]
git -C "$EXECUTOR" diff --quiet HEAD -- benchmark_tools
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
unset PYTHONPATH PYTHONHOME LD_PRELOAD LD_LIBRARY_PATH
cd "$ROOT"
exec /home/bizon/anaconda3/bin/python -B "$EXECUTOR/benchmark_tools/run_corrected_swiss_comparison.py" \
  --comparison "$ROOT/benchmark_tools/results/qfo_corrected_comparison_20260926_v7/manifest.json" \
  --comparison-sha256 042aa221d01554114a1b0b413ca3e2ff56dad5f8c06362c69a1bf49f2de642fc \
  --baseline "$ROOT/benchmark_tools/results/qfo_swiss_counts_20260917.json" \
  --protocol "$ROOT/benchmark_tools/results/QFO_SWISS_COMPARATOR_UNCERTAINTY_PROTOCOL_20260917.md" \
  --release-protocol "$ROOT/benchmark_tools/results/QFO_CORRECTED_RELEASE_PROTOCOL_20260917.md" \
  --output "$ROOT/benchmarks/work/qfo_recovered_swiss_uncertainty_${SLURM_JOB_ID}.json"
