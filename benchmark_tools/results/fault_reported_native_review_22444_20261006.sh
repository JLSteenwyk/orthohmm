#!/bin/bash
#SBATCH --job-name=ohmm_review08_fault_reported
#SBATCH --partition=gpu
#SBATCH --nodelist=bizon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=128G
#SBATCH --time=06:00:00
#SBATCH --no-requeue
set -euo pipefail
test "$#" -eq 2
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
cd "$ROOT"
exec env -u PYTHONHOME -u PYTHONPATH -u PYTHONUSERBASE \
    -u LD_PRELOAD -u LD_LIBRARY_PATH -u LD_AUDIT \
    PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 \
    PYTHONFAULTHANDLER=1 PYTHONUNBUFFERED=1 \
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    "$ROOT/benchmarks/work/native_factorial_review_py310_20261004/bin/python" \
    -B -X faulthandler -m benchmark_tools.run_fault_reported_native_review \
    --diagnostic-sha256 "$1" --worker-sha256 "$2"
