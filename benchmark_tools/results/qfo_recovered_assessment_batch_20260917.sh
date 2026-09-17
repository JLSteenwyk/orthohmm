#!/bin/bash
#SBATCH --job-name=qfo_recovered_assessment
#SBATCH --array=0-3%1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_recovered_assessment_%A_%a.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Frozen executor required}
COMMIT=${2:?Frozen executor revision required}
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --quiet HEAD -- benchmark_tools
cd "$ROOT"
exec /home/bizon/anaconda3/bin/python "$EXECUTOR/benchmark_tools/run_qfo_recovered_assessment.py" \
    --root "$ROOT" --index "$SLURM_ARRAY_TASK_ID" \
    --environment-sha256 e86545fd04cb644ed642cf2a31b4993fb2225ca4c05743cd03e661db129e82bc \
    --pairs-sha256 ce6f19cd005b886a91dc3aee63cb7b41e7f448d8cad5e65ff7cb10d92a636100
