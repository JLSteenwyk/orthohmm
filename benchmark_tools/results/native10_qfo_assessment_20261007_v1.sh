#!/bin/bash
#SBATCH --job-name=ohmm_native10_qfo_score
#SBATCH --partition=gpu
#SBATCH --nodelist=bizon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --no-requeue
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
cd "$ROOT"
if [[ "${SLURM_CPUS_PER_TASK:-}" != 8 || ! "${SLURM_JOB_ID:-}" =~ ^[1-9][0-9]*$ ]]; then
    printf '%s\n' 'Require a distinct scheduled eight-CPU assessment.' >&2
    exit 2
fi
exec env -u PYTHONHOME -u PYTHONPATH -u PYTHONUSERBASE \
    -u LD_PRELOAD -u LD_LIBRARY_PATH -u LD_AUDIT \
    PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 \
    PYTHONFAULTHANDLER=1 PYTHONUNBUFFERED=1 \
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    "$ROOT/benchmarks/work/native_factorial_review_py310_20261004/bin/python" \
    -B -X faulthandler -m benchmark_tools.run_allocated_native_factorial_qfo_assessment \
    --root "$ROOT" \
    --pairs "$ROOT/benchmarks/work/native10_qfo_pairs_allocated_20261007_v1/results.json" \
    --pairs-sha256 593679a7c35a6c5e2eb706846a6fcced706dd0d624bccb5c1608e6f83f7e7040 \
    --conversion-job 23977
