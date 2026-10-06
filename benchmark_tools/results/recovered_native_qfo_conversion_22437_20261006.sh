#!/bin/bash
#SBATCH --job-name=ohmm_qfo_recovered_pairs07
#SBATCH --partition=gpu
#SBATCH --nodelist=bizon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --no-requeue
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
cd "$ROOT"
exec env -u PYTHONHOME -u PYTHONPATH -u LD_PRELOAD -u LD_LIBRARY_PATH -u LD_AUDIT \
    PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 \
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    "$ROOT/benchmarks/work/native_factorial_review_py310_20261004/bin/python" -B \
    -m benchmark_tools.prepare_measurement_failed_native_qfo_pairs \
    --request "$ROOT/benchmarks/work/native_factorial_launch_20261004/request_07_receipt_amended.json" \
    --request-sha256 1078eb4bb39bd2f068e3faae3df1ca0131c2cf848c333aa0bc18424fe82213c7 \
    --scientific-recovery "$ROOT/benchmarks/work/native_factorial_measurement_failure_review_22437/review.json" \
    --scientific-recovery-sha256 650b69e0965e0a667359efaa43c1d814108f613547cbe66cf12964bd1e64e459 \
    --output-directory "$ROOT/benchmarks/work/measurement_failed_native_qfo_pairs_22437"
