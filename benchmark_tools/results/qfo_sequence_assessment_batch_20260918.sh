#!/bin/bash
#SBATCH --job-name=qfo_sequence_assess
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_sequence_assess_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Frozen executor required}
COMMIT=${2:?Exact revision required}
VARIANT=${3:?Variant required}
JOB=${4:?Successful conversion job required}
[[ $# == 4 && "$COMMIT" =~ ^[0-9a-f]{40}$ && "$JOB" =~ ^[0-9]+$ ]]
[[ "$VARIANT" == all_hits || "$VARIANT" == top100 ]]
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --exit-code HEAD -- benchmark_tools
MANIFEST="$ROOT/benchmarks/results/qfo_sequence_pairs_v1/$VARIANT/results.json"
read -r SHA _ < <(sha256sum "$MANIFEST")
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$ROOT"
exec /home/bizon/anaconda3/bin/python "$EXECUTOR/benchmark_tools/run_qfo_sequence_assessment.py" \
    --root "$ROOT" --variant "$VARIANT" --pairs-sha256 "$SHA" --conversion-job "$JOB"
