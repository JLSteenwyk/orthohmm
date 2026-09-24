#!/bin/bash
#SBATCH --job-name=qfo_recovery_score_admit
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_recovery_score_admit_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Frozen validator executor required}
COMMIT=${2:?Exact validator revision required}
SCORER=${3:?Frozen assessment executor required}
SCORER_COMMIT=${4:?Exact assessment revision required}
CONVERTER=${5:?Frozen conversion executor required}
CONVERTER_COMMIT=${6:?Exact conversion revision required}
JOB=${7:?Completed assessment job required}
CONVERSION_JOB=${8:?Completed conversion job required}
SHA=${9:?Pinned conversion manifest SHA256 required}
[[ "$JOB" =~ ^[0-9]+$ && "$CONVERSION_JOB" =~ ^[0-9]+$ && "$SHA" =~ ^[0-9a-f]{64}$ ]]
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --quiet HEAD -- benchmark_tools orthohmm
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$ROOT"
exec /home/bizon/anaconda3/bin/python "$EXECUTOR/benchmark_tools/admit_qfo_recovered_orthomcl_assessment.py" \
  --root "$ROOT" --job "$JOB" --conversion-job "$CONVERSION_JOB" --pairs-sha256 "$SHA" \
  --executor "$SCORER" --commit "$SCORER_COMMIT" --converter "$CONVERTER" --converter-commit "$CONVERTER_COMMIT" \
  --output "$ROOT/benchmarks/work/qfo_blast_recovery_score_admission_20260923.json"
