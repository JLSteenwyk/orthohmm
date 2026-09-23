#!/bin/bash
#SBATCH --job-name=qfo_corrected_fastoma_pairs
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_corrected_fastoma_pairs_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Pinned conversion executor required}
COMMIT=${2:?Exact executor revision required}
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --exit-code HEAD -- benchmark_tools orthohmm
ADMISSION="$ROOT/benchmarks/work/qfo_corrected_fastoma_admission_22054.json"
SHA=3cd4aac88a6e227493ab963b6a83a248adcde4780dc98d71ba3c401489805b8c
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
unset PYTHONNOUSERSITE PYTHONPATH PYTHONHOME LD_PRELOAD LD_LIBRARY_PATH
cd "$ROOT"
exec /home/bizon/anaconda3/bin/python -B "$EXECUTOR/benchmark_tools/prepare_qfo_corrected_fastoma_pairs.py" \
    --root "$ROOT" --admission "$ADMISSION" --admission-sha256 "$SHA" --admission-job 22054
