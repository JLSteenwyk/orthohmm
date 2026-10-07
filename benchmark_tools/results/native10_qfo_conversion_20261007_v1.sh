#!/bin/bash
#SBATCH --job-name=ohmm_native10_qfo_pairs
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
if [[ "${SLURM_CPUS_PER_TASK:-}" != 2 || ! "${SLURM_JOB_ID:-}" =~ ^[1-9][0-9]*$ ]]; then
    printf '%s\n' 'Require a distinct scheduled two-CPU conversion.' >&2
    exit 2
fi
exec env -u PYTHONHOME -u PYTHONPATH -u PYTHONUSERBASE \
    -u LD_PRELOAD -u LD_LIBRARY_PATH -u LD_AUDIT \
    PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 \
    PYTHONFAULTHANDLER=1 PYTHONUNBUFFERED=1 \
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    "$ROOT/benchmarks/work/native_factorial_review_py310_20261004/bin/python" \
    -B -X faulthandler -m benchmark_tools.prepare_allocated_native_factorial_qfo_pairs \
    --request "$ROOT/benchmarks/work/native_factorial_launch_20261004/request_10_allocated_v1.json" \
    --request-sha256 1355ae3b73a9c8496133bd25a6cbf136d149598ccab680672e81f86482338910 \
    --terminal-review "$ROOT/benchmarks/work/native10_review_finalization_20261007_v1/review/review.json" \
    --terminal-review-sha256 7eab212d9788deb031232e4641089ff96503c0e5727e1ec369067f93ffb518b5 \
    --output-directory "$ROOT/benchmarks/work/native10_qfo_pairs_allocated_20261007_v1"
