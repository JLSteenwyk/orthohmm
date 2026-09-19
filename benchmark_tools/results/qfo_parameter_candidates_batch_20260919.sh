#!/bin/bash
#SBATCH --job-name=qfo_parameter_candidates
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_parameter_candidates_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Frozen executor required}
COMMIT=${2:?Exact executor commit required}
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --quiet HEAD -- benchmark_tools orthohmm
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$ROOT"
exec /usr/bin/time -v -o "$ROOT/benchmarks/work/qfo_parameter_candidates_${SLURM_JOB_ID}.time" \
    /home/bizon/anaconda3/bin/python "$EXECUTOR/benchmark_tools/prepare_qfo_candidate_neighborhood.py" \
    --root "$ROOT" --output "$ROOT/benchmarks/results/qfo_parameter_candidates_v1"
