#!/bin/bash
#SBATCH --job-name=qfo_corrected_of_assess
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --array=0-1%1
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_corrected_of_assess_%A_%a.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Pinned executor required}
COMMIT=${2:?Exact revision required}
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --exit-code HEAD -- benchmark_tools
case "$SLURM_ARRAY_TASK_ID" in
    0) METHOD=orthofinder_full ;;
    1) METHOD=orthofinder_sequence_only ;;
    *) exit 2 ;;
esac
PAIRS="$ROOT/benchmarks/results/qfo_corrected_comparator_pairs_v1/$METHOD/results.json"
read -r SHA _ < <(sha256sum "$PAIRS")
JOB=$(jq -er '.job_id | strings' "$PAIRS")
[[ "$JOB" =~ ^[0-9]+$ ]]
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$ROOT"
exec /home/bizon/anaconda3/bin/python "$EXECUTOR/benchmark_tools/run_qfo_corrected_comparator_assessment.py" \
    --root "$ROOT" --method "$METHOD" --pairs-sha256 "$SHA" --conversion-job "$JOB"
