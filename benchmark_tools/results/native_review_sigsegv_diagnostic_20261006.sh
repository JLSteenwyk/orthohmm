#!/bin/bash
#SBATCH --job-name=ohmm_review08_diagnostic
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
DEST="$ROOT/benchmarks/work/native_review_sigsegv_diagnostic_20261006_v1"
cd "$ROOT"
test ! -e "$DEST"
test ! -L "$DEST"
mkdir "$DEST"
exec /usr/bin/time -v -o "$DEST/time.txt" env \
    -u PYTHONHOME -u PYTHONPATH -u PYTHONUSERBASE \
    -u LD_PRELOAD -u LD_LIBRARY_PATH -u LD_AUDIT \
    PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 \
    PYTHONFAULTHANDLER=1 PYTHONUNBUFFERED=1 \
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    "$ROOT/benchmarks/work/native_factorial_review_py310_20261004/bin/python" \
    -B -X faulthandler -m benchmark_tools.validate_native_factorial_outputs \
    --request "$ROOT/benchmarks/work/native_factorial_launch_20261004/request_08_receipt_amended.json" \
    --request-sha256 35f8ef9b1d7abc1f574b8e7e3c55bc91c68c9f8f74247a9a816754a41a93fb6a \
    --output "$DEST/outputs.json"
