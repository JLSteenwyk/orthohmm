#!/bin/bash
#SBATCH --job-name=qfo_recovery_orthomcl_assess
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_recovery_assess_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Frozen assessment executor required}
COMMIT=${2:?Exact assessment revision required}
CONVERTER=${3:?Frozen conversion executor required}
CONVERTER_COMMIT=${4:?Exact conversion revision required}
JOB=${5:?Completed conversion job required}
SHA=${6:?Pinned conversion manifest SHA256 required}
[[ "$JOB" =~ ^[0-9]+$ && "$SHA" =~ ^[0-9a-f]{64}$ ]]
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --quiet HEAD -- benchmark_tools orthohmm
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$ROOT"
exec /home/bizon/anaconda3/bin/python "$EXECUTOR/benchmark_tools/run_qfo_recovered_orthomcl_assessment.py" \
  --root "$ROOT" --pairs-sha256 "$SHA" --conversion-job "$JOB" \
  --executor "$CONVERTER" --commit "$CONVERTER_COMMIT"
