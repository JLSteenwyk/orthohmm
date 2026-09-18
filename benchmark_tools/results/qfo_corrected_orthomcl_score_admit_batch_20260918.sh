#!/bin/bash
#SBATCH --job-name=qfo_corrected_orthomcl_score_admit
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_corrected_orthomcl_score_admit_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Pinned executor required}
COMMIT=${2:?Exact revision required}
CONVERSION_JOB=${3:?Conversion job required}
JOB=${4:?Scoring job required}
[[ "$JOB" =~ ^[0-9]+$ && "$CONVERSION_JOB" =~ ^[0-9]+$ ]]
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --exit-code HEAD -- benchmark_tools
PAIRS="$ROOT/benchmarks/results/qfo_corrected_comparator_pairs_v1/orthomcl/results.json"
read -r SHA _ < <(sha256sum "$PAIRS")
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$ROOT"
exec /home/bizon/anaconda3/bin/python "$EXECUTOR/benchmark_tools/admit_qfo_corrected_comparator_assessment.py" \
    --root "$ROOT" --method orthomcl --job "$JOB" --conversion-job "$CONVERSION_JOB" \
    --pairs-sha256 "$SHA" --output "$ROOT/benchmarks/work/qfo_corrected_orthomcl_assessment_admission_20260918.json"
