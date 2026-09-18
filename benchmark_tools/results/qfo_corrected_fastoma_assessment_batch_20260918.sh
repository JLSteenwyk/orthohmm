#!/bin/bash
#SBATCH --job-name=qfo_corrected_fastoma_assess
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_corrected_fastoma_assess_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Pinned assessment executor required}
COMMIT=${2:?Exact executor revision required}
CONVERSION_JOB=${3:?Pair conversion job required}
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --exit-code HEAD -- benchmark_tools orthohmm
PAIRS="$ROOT/benchmarks/results/qfo_corrected_comparator_pairs_v1/fastoma/results.json"
read -r SHA _ < <(sha256sum "$PAIRS")
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$ROOT"
exec /home/bizon/anaconda3/bin/python "$EXECUTOR/benchmark_tools/run_qfo_corrected_comparator_assessment.py" \
    --root "$ROOT" --method fastoma --pairs-sha256 "$SHA" --conversion-job "$CONVERSION_JOB"
