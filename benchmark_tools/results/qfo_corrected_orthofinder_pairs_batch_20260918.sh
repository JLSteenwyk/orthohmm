#!/bin/bash
#SBATCH --job-name=qfo_corrected_of_pairs
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --array=0-1%1
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_corrected_of_pairs_%A_%a.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Pinned executor required}
COMMIT=${2:?Exact revision required}
ADMISSION_JOB=${3:?Native admission job required}
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --exit-code HEAD -- benchmark_tools orthohmm
case "$SLURM_ARRAY_TASK_ID" in
    0) METHOD=orthofinder_full ;;
    1) METHOD=orthofinder_sequence_only ;;
    *) exit 2 ;;
esac
ADMISSION="$ROOT/benchmarks/work/qfo_corrected_orthofinder_admission_20260918.json"
read -r SHA _ < <(sha256sum "$ADMISSION")
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$ROOT"
exec /home/bizon/anaconda3/bin/python "$EXECUTOR/benchmark_tools/prepare_qfo_corrected_orthofinder_pairs.py" \
    --root "$ROOT" --method "$METHOD" --admission "$ADMISSION" \
    --admission-sha256 "$SHA" --admission-job "$ADMISSION_JOB"
