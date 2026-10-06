#!/bin/bash
#SBATCH --job-name=ohmm_recover07
#SBATCH --partition=gpu
#SBATCH --nodelist=bizon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/native_factorial_launch_20261004/recovery07_%j.out
#SBATCH --error=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/native_factorial_launch_20261004/recovery07_%j.err
set -euo pipefail
cd /mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
exec env -u PYTHONHOME -u PYTHONPATH -u LD_PRELOAD -u LD_LIBRARY_PATH -u LD_AUDIT \
  PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  benchmarks/work/native_factorial_review_py310_20261004/bin/python -B \
  -m benchmark_tools.review_native_factorial_measurement_failure \
  --request benchmarks/work/native_factorial_launch_20261004/request_07_receipt_amended.json \
  --request-sha256 1078eb4bb39bd2f068e3faae3df1ca0131c2cf848c333aa0bc18424fe82213c7 \
  --audit benchmark_tools/results/native_factorial_cadence_failure_22437_20261006.json \
  --audit-sha256 78a13590e1d33e89c91a80c88ad41e561d1475f2c1ce1cbb1663df5a99949e99 \
  --failed-review benchmarks/work/native_factorial_terminal_review_22437/failure.json \
  --failed-review-sha256 2a11c748611dd01b3c29b86bf072be9bef5f864fb30e8d6dcd91fba11de81fbb \
  --output-directory benchmarks/work/native_factorial_measurement_failure_review_22437
