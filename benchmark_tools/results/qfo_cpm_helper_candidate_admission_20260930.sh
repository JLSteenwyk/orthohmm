#!/bin/bash
#SBATCH --job-name=qfo_cpm_helper_candidate_admission
#SBATCH --nodelist=bizon
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_cpm_helper_candidate_admission_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Frozen executor required}
COMMIT=${2:?Exact executor commit required}
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --quiet HEAD -- benchmark_tools orthohmm
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
unset PYTHONPATH PYTHONHOME LD_PRELOAD LD_LIBRARY_PATH
cd "$ROOT"
OUTPUT="$ROOT/benchmarks/work/qfo_cpm_helper_candidate_admission_20260930"
mkdir "$OUTPUT"
exec /usr/bin/time -v -o "$ROOT/benchmarks/work/qfo_cpm_helper_candidate_admission_${SLURM_JOB_ID}.time.txt" \
  /home/bizon/anaconda3/bin/python -B "$EXECUTOR/benchmark_tools/admit_helper_cpm_candidates.py" \
  --root "$ROOT" --output "$OUTPUT/status.json" \
  --preparation-sha256 2e62321bf7a6d410f8ddafff783aaea979f1fec753fd6415b2b439ee7e2f54e6 \
  --protocol-sha256 96f3cd4baaea3cfaf4fb949d2f683888bac41b8b2e7696c37cc85f747cdd5955
