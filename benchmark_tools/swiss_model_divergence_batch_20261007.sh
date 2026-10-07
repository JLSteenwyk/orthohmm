#!/bin/bash
#SBATCH --job-name=ohmm-swiss-divergence
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --no-requeue
set -euo pipefail
ROOT="${SLURM_SUBMIT_DIR:?}"
exec env -u PYTHONPATH -u PYTHONHOME -u PYTHONUSERBASE \
    -u LD_PRELOAD -u LD_LIBRARY_PATH -u LD_AUDIT \
    PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 \
    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    "$ROOT/benchmarks/work/native_factorial_review_py310_20261004/bin/python" \
    -B -m benchmark_tools.prepare_swiss_model_divergence \
    --repo "$ROOT" --source-commit "${1:?}" --output "${2:?}" \
    --job-id "${SLURM_JOB_ID:?}"
