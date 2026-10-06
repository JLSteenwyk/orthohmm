#!/bin/bash
#SBATCH --job-name=ohmm_qfo_native_assess08
#SBATCH --partition=gpu
#SBATCH --nodelist=bizon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --no-requeue
set -euo pipefail
cd /mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
exec env -u PYTHONHOME -u PYTHONPATH -u LD_PRELOAD -u LD_LIBRARY_PATH -u LD_AUDIT \
  PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  benchmarks/work/native_factorial_review_py310_20261004/bin/python -B \
  -m benchmark_tools.run_conversion_gated_native_qfo_assessment \
  --conversion-submission benchmark_tools/results/review_gated_native_qfo_conversion_submission_22450.json \
  --conversion-submission-sha256 9fab579abb0de79ef9cb163a4afb3aad76776021fe1eb86c9fe8f3b828af3e4f \
  --worker-sha256 797902b45f3decc1b3c20ce03fc4c9084bce8092efa9f0e66b4a13be3f80fbde
