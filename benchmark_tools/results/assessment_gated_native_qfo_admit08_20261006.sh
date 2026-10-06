#!/bin/bash
#SBATCH --job-name=ohmm_qfo_native_admit08
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
  -m benchmark_tools.run_assessment_gated_native_qfo_admission \
  --assessment-submission benchmark_tools/results/conversion_gated_native_qfo_assessment_submission_22451.json \
  --assessment-submission-sha256 ce8cd608e3a0d2cb62fcda052c5c14bef6b81d796ff65967f0df98ea11c517b6 \
  --worker-sha256 ba6c0d191e9c5582a0594c707cc1785c6d66b2433dd0803e16fee176abd690a6
