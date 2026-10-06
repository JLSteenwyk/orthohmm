#!/bin/bash
#SBATCH --job-name=ohmm_qfo_fault_review_pairs08
#SBATCH --partition=gpu
#SBATCH --nodelist=bizon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --no-requeue
set -euo pipefail
if [[ $# -ne 1 || ! ${1:-} =~ ^[0-9a-f]{64}$ ]]; then
    printf '%s\n' 'Require the independently verified terminal-review SHA256.' >&2
    exit 2
fi
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
cd "$ROOT"
driver_digest=$(sha256sum "$ROOT/benchmark_tools/prepare_native_factorial_qfo_pairs.py")
test "${driver_digest%% *}" = 472f26a140cf89aed36f8212a31fef7d8f5061204dac2e0ce80c3c7e6746630e
test "${SLURM_CPUS_PER_TASK:-}" = 2
test "${SLURM_JOB_ID:-}" != 23017
test "${SLURM_JOB_ID:-}" != 22444
test -n "${SLURM_JOB_ID:-}"
exec env -u PYTHONHOME -u PYTHONPATH -u PYTHONUSERBASE \
    -u LD_PRELOAD -u LD_LIBRARY_PATH -u LD_AUDIT \
    PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 \
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    "$ROOT/benchmarks/work/native_factorial_review_py310_20261004/bin/python" -B \
    -m benchmark_tools.prepare_native_factorial_qfo_pairs \
    --request "$ROOT/benchmarks/work/native_factorial_launch_20261004/request_08_receipt_amended.json" \
    --request-sha256 35f8ef9b1d7abc1f574b8e7e3c55bc91c68c9f8f74247a9a816754a41a93fb6a \
    --terminal-review "$ROOT/benchmarks/work/native_factorial_terminal_review_22444_fault_reported_v1/review.json" \
    --terminal-review-sha256 "$1" \
    --output-directory "$ROOT/benchmarks/work/native_factorial_qfo_pairs_22444_fault_reported_v1"
