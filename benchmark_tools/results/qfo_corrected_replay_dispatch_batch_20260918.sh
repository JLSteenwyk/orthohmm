#!/bin/bash
#SBATCH --job-name=qfo_corrected_checked_replay
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=32
#SBATCH --mem=192G
#SBATCH --time=24:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_corrected_checked_replay_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Frozen dispatch executor required}
COMMIT=${2:?Exact executor revision required}
PREPARATION_JOB=${3:?Plan preparation job required}
ADMISSION_JOB=${4:?HMM admission job required}
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --exit-code HEAD -- benchmark_tools orthohmm
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$ROOT"
exec /home/bizon/anaconda3/bin/python "$EXECUTOR/benchmark_tools/launch_qfo_corrected_replay.py" \
    --root "$ROOT" --preparation-job "$PREPARATION_JOB" --admission-job "$ADMISSION_JOB"
