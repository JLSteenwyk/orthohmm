#!/bin/bash
#SBATCH --job-name=qfo_recovery_pairs
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_recovery_pairs_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Frozen conversion executor required}
COMMIT=${2:?Exact conversion revision required}
ADMITTER=${3:?Frozen native admission executor required}
ADMITTER_COMMIT=${4:?Exact native admission revision required}
JOB=${5:?Completed native admission job required}
SHA=${6:?Pinned native admission report SHA256 required}
[[ "$JOB" =~ ^[0-9]+$ && "$SHA" =~ ^[0-9a-f]{64}$ ]]
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --quiet HEAD -- benchmark_tools orthohmm
CACHE="/tmp/orthohmm-recovery-pairs-cache-${SLURM_JOB_ID:?}"
[[ ! -e "$CACHE" && ! -L "$CACHE" ]]
cd "$ROOT"
exec env -i HOME="$HOME" USER=bizon LOGNAME=bizon PATH=/usr/bin:/bin LANG=C LC_ALL=C \
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  SLURM_JOB_ID="$SLURM_JOB_ID" SLURM_CPUS_PER_TASK="$SLURM_CPUS_PER_TASK" \
  SLURM_MEM_PER_NODE="$SLURM_MEM_PER_NODE" \
  "$ROOT/benchmarks/work/orthomcl_python_env_20260918/bin/python" -I -B -X "pycache_prefix=$CACHE" \
  "$EXECUTOR/benchmark_tools/prepare_recovered_orthomcl_pairs.py" --root "$ROOT" \
  --admission "$ROOT/benchmarks/results/qfo_blast_recovery_native_admission_v1/report.json" \
  --admission-sha256 "$SHA" --job "$JOB" --executor "$ADMITTER" --commit "$ADMITTER_COMMIT"
