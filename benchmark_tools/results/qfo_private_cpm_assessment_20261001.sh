#!/bin/bash
#SBATCH --job-name=qfo_private_cpm_assessment
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_private_cpm_assessment_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Frozen assessment executor required}
COMMIT=${2:?Exact assessment revision required}
PROTOCOL_SHA=${3:?Reviewed assessment protocol required}
CONVERSION_JOB=${4:?Completed pair conversion job required}
CONVERSION_SHA=${5:?Reviewed pair conversion report required}
SUBMISSION_SHA=${6:?Reviewed pair conversion submission required}
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --quiet HEAD -- benchmark_tools orthohmm
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
unset PYTHONNOUSERSITE PYTHONPATH PYTHONHOME LD_PRELOAD LD_LIBRARY_PATH LD_AUDIT
cd "$ROOT"
exec /home/bizon/anaconda3/bin/python -B "$EXECUTOR/benchmark_tools/run_private_helper_cpm_assessment.py" \
  --root "$ROOT" --protocol-sha256 "$PROTOCOL_SHA" \
  --conversion-job "$CONVERSION_JOB" --conversion-sha256 "$CONVERSION_SHA" \
  --conversion-submission-sha256 "$SUBMISSION_SHA"
