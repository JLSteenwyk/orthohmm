#!/bin/bash
#SBATCH --job-name=qfo_corrected_bpo_prepare
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_corrected_bpo_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Pinned executor required}
COMMIT=${2:?Exact executor revision required}
ADMISSION_JOB=${3:?BLAST admission job required}
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --exit-code HEAD -- benchmark_tools orthohmm
ADMISSION="$ROOT/benchmarks/work/qfo_corrected_blast_admission_20260918/report.json"
read -r SHA _ < <(sha256sum "$ADMISSION")
CACHE="/tmp/orthohmm-bpo-cache-${SLURM_JOB_ID:?}"
[[ ! -e "$CACHE" && ! -L "$CACHE" ]]
cd "$ROOT"
exec env -i HOME="$HOME" USER=bizon LOGNAME=bizon PATH=/usr/bin:/bin LANG=C LC_ALL=C \
    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    SLURM_JOB_ID="$SLURM_JOB_ID" SLURM_CPUS_PER_TASK="$SLURM_CPUS_PER_TASK" \
    SLURM_MEM_PER_NODE="$SLURM_MEM_PER_NODE" \
    "$ROOT/benchmarks/work/orthomcl_python_env_20260918/bin/python" -I -B -X "pycache_prefix=$CACHE" \
    "$EXECUTOR/benchmark_tools/prepare_qfo_corrected_bpo.py" \
    --root "$ROOT" --admission "$ADMISSION" --admission-sha256 "$SHA" --admission-job "$ADMISSION_JOB"
