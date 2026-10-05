#!/bin/bash
#SBATCH --job-name=orthohmm_native_review
#SBATCH --nodelist=bizon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --no-requeue
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
unset PYTHONHOME PYTHONPATH LD_PRELOAD LD_LIBRARY_PATH LD_AUDIT
PYTHON="$ROOT/benchmarks/work/native_factorial_review_py310_20261004/bin/python"
read -r ACTUAL _ < <(sha256sum -- "$PYTHON")
[[ "$ACTUAL" == 8b1cd756be711ef53f35cb6c954472fdfc52094c4619a01f553f64354587388b ]]
cd "$ROOT"
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
exec "$PYTHON" -B benchmark_tools/finish_native_orthobench_attempt.py "$@"
