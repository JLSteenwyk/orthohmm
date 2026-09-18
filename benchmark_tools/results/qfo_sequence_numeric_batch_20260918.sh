#!/bin/bash
#SBATCH --job-name=qfo_corrected_sequence_numeric
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=2
#SBATCH --mem=192G
#SBATCH --time=7-00:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_sequence_numeric_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Frozen converter required}
COMMIT=${2:?Exact revision required}
JOB=${3:?Search admission job required}
[[ $# == 3 && "$COMMIT" =~ ^[0-9a-f]{40}$ && "$JOB" =~ ^[0-9]+$ ]]
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --exit-code HEAD -- benchmark_tools orthohmm
ADMISSION="$ROOT/benchmarks/work/qfo_sequence_search_admission_20260918.json"
read -r SHA _ < <(sha256sum "$ADMISSION")
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$ROOT"
exec /home/bizon/anaconda3/bin/python "$EXECUTOR/benchmark_tools/convert_qfo_sequence_search_control.py" \
    --root "$ROOT" --admission "$ADMISSION" --admission-sha256 "$SHA" --admission-job "$JOB" \
    --output "$ROOT/benchmarks/results/qfo_sequence_numeric_v1"
