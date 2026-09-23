#!/bin/bash
#SBATCH --job-name=qfo_corrected_fastoma_admit
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_corrected_fastoma_admit_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Pinned admission executor required}
COMMIT=${2:?Exact executor revision required}
INFERENCE_JOB=${3:?Native inference job required}
STAGE_JOB=${4:?Input staging job required}
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --exit-code HEAD -- benchmark_tools orthohmm
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
unset PYTHONNOUSERSITE PYTHONPATH PYTHONHOME LD_PRELOAD LD_LIBRARY_PATH
cd "$ROOT"
exec /home/bizon/anaconda3/bin/python -B "$EXECUTOR/benchmark_tools/admit_qfo_corrected_fastoma.py" \
    --root "$ROOT" --job "$INFERENCE_JOB" --stage-job "$STAGE_JOB" \
    --output "$ROOT/benchmarks/work/qfo_corrected_fastoma_admission_${SLURM_JOB_ID}.json"
