#!/bin/bash
#SBATCH --job-name=qfo_private_cpm_pair_conversion
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_private_cpm_pair_conversion_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Frozen conversion executor required}
COMMIT=${2:?Exact conversion revision required}
PROTOCOL_SHA=${3:?Reviewed conversion protocol required}
ADMISSION_JOB=${4:?Completed native admission job required}
ADMISSION_SHA=${5:?Reviewed native admission report required}
SUBMISSION_SHA=${6:?Reviewed native admission submission required}
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --quiet HEAD -- benchmark_tools orthohmm
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
unset PYTHONNOUSERSITE PYTHONPATH PYTHONHOME LD_PRELOAD LD_LIBRARY_PATH LD_AUDIT
cd "$ROOT"
exec /home/bizon/anaconda3/bin/python -B "$EXECUTOR/benchmark_tools/prepare_private_helper_cpm_pairs.py" \
  --root "$ROOT" --protocol-sha256 "$PROTOCOL_SHA" \
  --admission-job "$ADMISSION_JOB" --admission-sha256 "$ADMISSION_SHA" \
  --admission-submission-sha256 "$SUBMISSION_SHA"
