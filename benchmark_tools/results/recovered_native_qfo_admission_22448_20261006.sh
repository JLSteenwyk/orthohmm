#!/bin/bash
#SBATCH --job-name=ohmm_qfo_recovered_admit07
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
    -m benchmark_tools.admit_measurement_failed_native_qfo_assessment \
    --root "$ROOT" \
    --pairs "$ROOT/benchmarks/work/measurement_failed_native_qfo_pairs_22437_snapshot_restored/results.json" \
    --pairs-sha256 5b29ded11bd81b2b077b128a373ab2d1da037c5886b3a005ff8970bc09ca169c \
    --conversion-job 22447 --assessment-job 22448 \
    --output-directory "$ROOT/benchmarks/results/measurement_failed_native_qfo_admission_v1/p0_c0_r1"
