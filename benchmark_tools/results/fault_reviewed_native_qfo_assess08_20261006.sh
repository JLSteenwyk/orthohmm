#!/bin/bash
#SBATCH --job-name=ohmm_qfo_fault_review_assess08
#SBATCH --partition=gpu
#SBATCH --nodelist=bizon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --no-requeue
set -euo pipefail
if [[ $# -ne 2 || ! ${1:-} =~ ^[0-9a-f]{64}$ || ! ${2:-} =~ ^[1-9][0-9]*$ ]]; then
    printf '%s\n' 'Require the verified conversion SHA256 and completed conversion job ID.' >&2
    exit 2
fi
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
cd "$ROOT"
driver_digest=$(sha256sum "$ROOT/benchmark_tools/run_native_factorial_qfo_assessment.py")
test "${driver_digest%% *}" = ba051efa4b2dd9f72a4464f56a891b4742fa3a1c6999f8b95a65d0805c77dafe
test "${SLURM_CPUS_PER_TASK:-}" = 8
test "${SLURM_JOB_ID:-}" != "$2"
test "${SLURM_JOB_ID:-}" != 23017
test "${SLURM_JOB_ID:-}" != 22444
test -n "${SLURM_JOB_ID:-}"
exec env -u PYTHONHOME -u PYTHONPATH -u PYTHONUSERBASE \
    -u LD_PRELOAD -u LD_LIBRARY_PATH -u LD_AUDIT \
    PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 \
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    "$ROOT/benchmarks/work/native_factorial_review_py310_20261004/bin/python" -B \
    -m benchmark_tools.run_native_factorial_qfo_assessment \
    --root "$ROOT" \
    --pairs "$ROOT/benchmarks/work/native_factorial_qfo_pairs_22444_fault_reported_v1/results.json" \
    --pairs-sha256 "$1" --conversion-job "$2"
