#!/bin/bash
#SBATCH --job-name=ohmm_qfo_native_pairs08
#SBATCH --partition=gpu
#SBATCH --nodelist=bizon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --no-requeue
set -euo pipefail
cd /mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
exec env -u PYTHONHOME -u PYTHONPATH -u LD_PRELOAD -u LD_LIBRARY_PATH -u LD_AUDIT \
  PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  benchmarks/work/native_factorial_review_py310_20261004/bin/python -B \
  -m benchmark_tools.run_review_gated_native_qfo_pairs \
  --review-submission benchmark_tools/results/native_factorial_review_submission_22445.json \
  --review-submission-sha256 97cac14882db9ad15ae597e569a77ae8f5f5a6b748915948182d2cfae699f62b \
  --worker-sha256 bbb1681c6d68da076fac6d51815404282ce635fe1273eb05c5e6b040cb17e0ee
