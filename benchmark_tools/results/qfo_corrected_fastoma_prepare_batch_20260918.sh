#!/bin/bash
#SBATCH --job-name=qfo_corrected_fastoma_prepare
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_corrected_fastoma_prepare_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Pinned preparation executor required}
COMMIT=${2:?Exact executor revision required}
ADMISSION_JOB=${3:?Native OrthoFinder admission job required}
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --exit-code HEAD -- benchmark_tools orthohmm
ADMISSION="$ROOT/benchmarks/work/qfo_corrected_orthofinder_admission_20260918.json"
read -r SHA _ < <(sha256sum "$ADMISSION")
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$ROOT"
exec /home/bizon/anaconda3/bin/python "$EXECUTOR/benchmark_tools/prepare_qfo_corrected_fastoma.py" \
    --root "$ROOT" --admission "$ADMISSION" --admission-sha256 "$SHA" \
    --admission-job "$ADMISSION_JOB" \
    --directory "$ROOT/benchmarks/work/qfo_corrected_fastoma_inputs_20260918" \
    --output "$ROOT/benchmarks/work/qfo_corrected_fastoma_staging_20260918.json"
