#!/bin/bash
#SBATCH --job-name=qfo_recovery_native_admit
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_recovery_native_admit_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Frozen validator executor required}
COMMIT=${2:?Exact validator revision required}
NATIVE=${3:?Frozen native executor required}
NATIVE_COMMIT=${4:?Exact native revision required}
BPO=${5:?Frozen BPO admission executor required}
BPO_COMMIT=${6:?Exact BPO admission revision required}
JOB=${7:?Completed native job required}
[[ "$JOB" =~ ^[0-9]+$ ]]
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --quiet HEAD -- benchmark_tools orthohmm
CACHE="/tmp/orthohmm-recovery-native-admit-cache-${SLURM_JOB_ID:?}"
[[ ! -e "$CACHE" && ! -L "$CACHE" ]]
cd "$ROOT"
exec env -i HOME="$HOME" USER=bizon LOGNAME=bizon PATH=/usr/bin:/bin LANG=C LC_ALL=C \
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  SLURM_JOB_ID="$SLURM_JOB_ID" SLURM_CPUS_PER_TASK="$SLURM_CPUS_PER_TASK" \
  SLURM_MEM_PER_NODE="$SLURM_MEM_PER_NODE" \
  "$ROOT/benchmarks/work/orthomcl_python_env_20260918/bin/python" -I -B -X "pycache_prefix=$CACHE" \
  "$EXECUTOR/benchmark_tools/admit_recovered_orthomcl.py" --root "$ROOT" \
  --job "$JOB" --executor "$NATIVE" --commit "$NATIVE_COMMIT" \
  --bpo-executor "$BPO" --bpo-commit "$BPO_COMMIT" \
  --output "$ROOT/benchmarks/results/qfo_blast_recovery_native_admission_v1"
