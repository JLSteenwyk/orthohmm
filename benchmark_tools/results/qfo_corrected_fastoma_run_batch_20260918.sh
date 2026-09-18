#!/bin/bash
#SBATCH --job-name=qfo_corrected_fastoma
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=180
#SBATCH --mem=720G
#SBATCH --time=7-00:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_corrected_fastoma_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Pinned inference executor required}
COMMIT=${2:?Exact executor revision required}
STAGE_JOB=${3:?Successful staging job required}
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --exit-code HEAD -- benchmark_tools orthohmm
STAGE="$ROOT/benchmarks/work/qfo_corrected_fastoma_staging_20260918.json"
read -r SHA _ < <(sha256sum "$STAGE")
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$ROOT"
exec /home/bizon/anaconda3/bin/python "$EXECUTOR/benchmark_tools/run_qfo_corrected_fastoma.py" \
    --root "$ROOT" --stage "$STAGE" --stage-sha256 "$SHA" --stage-job "$STAGE_JOB"
