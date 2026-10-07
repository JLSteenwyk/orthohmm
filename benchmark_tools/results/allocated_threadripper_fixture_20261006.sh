#!/bin/bash
#SBATCH --job-name=ohmm_allocated_core_fixture
#SBATCH --partition=gpu
#SBATCH --nodelist=bizon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --time=00:20:00
#SBATCH --no-requeue
set -euo pipefail
if [[ $# -ne 2 || ! ${1:-} =~ ^[0-9a-f]{64}$ || ! ${2:-} =~ ^[0-9a-f]{64}$ ]]; then
    printf '%s\n' 'Require prepared fixture and driver SHA256 values.' >&2
    exit 2
fi
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
cd "$ROOT"
actual=$(sha256sum benchmark_tools/run_allocated_threadripper_fixture.py)
test "${actual%% *}" = "$2"
exec env -u PYTHONHOME -u PYTHONPATH -u PYTHONUSERBASE \
    -u LD_PRELOAD -u LD_LIBRARY_PATH -u LD_AUDIT \
    PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 \
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    "$ROOT/benchmarks/work/native_factorial_review_py310_20261004/bin/python" -B \
    -m benchmark_tools.run_allocated_threadripper_fixture \
    --prepared "$ROOT/benchmark_tools/results/allocated_threadripper_fixture_prepared_20261006.json" \
    --prepared-sha256 "$1"
