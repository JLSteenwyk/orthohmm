#!/bin/bash
#SBATCH --job-name=qfo_recovery_bpo_admit
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_recovery_bpo_admit_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Frozen admission executor required}
COMMIT=${2:?Exact admission revision required}
PREPARER=${3:?Frozen preparation executor required}
PREPARER_COMMIT=${4:?Exact preparation revision required}
JOB=${5:?Completed preparation job required}
[[ "$JOB" =~ ^[0-9]+$ ]]
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --quiet HEAD -- benchmark_tools orthohmm
CACHE="/tmp/orthohmm-recovery-bpo-admission-cache-${SLURM_JOB_ID:?}"
[[ ! -e "$CACHE" && ! -L "$CACHE" ]]
cd "$ROOT"
exec env -i HOME="$HOME" USER=bizon LOGNAME=bizon PATH=/usr/bin:/bin LANG=C LC_ALL=C \
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  "$ROOT/benchmarks/work/orthomcl_python_env_20260918/bin/python" -I -B -X "pycache_prefix=$CACHE" \
  "$EXECUTOR/benchmark_tools/admit_blast_recovery_bpo.py" --root "$ROOT" \
  --job "$JOB" --executor "$PREPARER" --commit "$PREPARER_COMMIT" \
  --output "$ROOT/benchmarks/results/qfo_blast_recovery_bpo_admission_v1"
