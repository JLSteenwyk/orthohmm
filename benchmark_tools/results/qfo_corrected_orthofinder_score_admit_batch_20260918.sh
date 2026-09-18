#!/bin/bash
#SBATCH --job-name=qfo_corrected_of_score_admit
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --array=0-1%1
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_corrected_of_score_admit_%A_%a.log
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
ASSESSMENT="$ROOT/benchmarks/results/qfo_corrected_assessment_v1/$METHOD/results.json"
read -r SHA _ < <(sha256sum "$PAIRS")
CONVERSION_JOB=$(jq -er '.job_id | strings' "$PAIRS")
JOB=$(jq -er '.job_id | strings' "$ASSESSMENT")
[[ "$JOB" =~ ^[0-9]+$ && "$CONVERSION_JOB" =~ ^[0-9]+$ ]]
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$ROOT"
exec /home/bizon/anaconda3/bin/python "$EXECUTOR/benchmark_tools/admit_qfo_corrected_comparator_assessment.py" \
    --root "$ROOT" --method "$METHOD" --job "$JOB" --conversion-job "$CONVERSION_JOB" \
    --pairs-sha256 "$SHA" --output "$ROOT/benchmarks/work/qfo_corrected_${METHOD}_assessment_admission_20260918.json"
