#!/bin/bash
#SBATCH --job-name=qfo_corrected_native_pairs
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_corrected_native_pairs_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Frozen converter executor required}
COMMIT=${2:?Exact revision required}
INDEX=${3:?R-on factorial index required}
CANDIDATE_JOB=${4:?Candidate admission job required}
[[ $# == 4 && "$COMMIT" =~ ^[0-9a-f]{40}$ && "$INDEX" =~ ^[1357]$ && "$CANDIDATE_JOB" =~ ^[0-9]+$ ]]
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --quiet HEAD -- benchmark_tools orthohmm qfo_benchmark/og_to_pairwise.py
CANDIDATE="$ROOT/benchmarks/work/qfo_corrected_candidate_admission_20260918.json"
NATIVE="$ROOT/benchmarks/work/qfo_corrected_factorial_native_admission_$((INDEX / 2))_20260918.json"
read -r CANDIDATE_SHA _ < <(sha256sum "$CANDIDATE")
read -r NATIVE_SHA _ < <(sha256sum "$NATIVE")
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$ROOT"
exec /home/bizon/anaconda3/bin/python "$EXECUTOR/benchmark_tools/prepare_qfo_corrected_native_pairs.py" \
    --root "$ROOT" --index "$INDEX" --candidate-admission "$CANDIDATE" \
    --candidate-sha256 "$CANDIDATE_SHA" --candidate-job "$CANDIDATE_JOB" \
    --native-admission "$NATIVE" --native-sha256 "$NATIVE_SHA"
