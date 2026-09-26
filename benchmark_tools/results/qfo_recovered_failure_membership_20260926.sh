#!/bin/bash
#SBATCH --job-name=qfo_recovered_failure_membership
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_recovered_failure_membership_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR="$ROOT/benchmarks/work/publication_native_failure_membership_v1_20260925"
COMMIT=9b9ab8dc7f56b7dfc29e2eec2f94cc919b06a018
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
[[ -z $(git -C "$EXECUTOR" status --porcelain) ]]
OUTPUT="$ROOT/benchmarks/work/qfo_recovered_failure_membership_20260926.json"
[[ ! -e "$OUTPUT" && ! -L "$OUTPUT" ]]
cd "$ROOT"
exec env -i HOME="$HOME" USER=bizon LOGNAME=bizon PATH=/usr/bin:/bin LANG=C LC_ALL=C \
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  "$ROOT/benchmarks/work/orthomcl_python_env_20260918/bin/python" -I -B \
  "$EXECUTOR/benchmark_tools/audit_orthomcl_failure_membership.py" \
  --search "$ROOT/benchmarks/results/qfo_blast_native_representation_admission_v1/report.json" \
  --search-sha256 b034aef886b4a68a5915f5796a2345fa10bd3668116d8de0b927d6f5bf00e4ff \
  --native "$ROOT/benchmarks/results/qfo_blast_recovery_native_admission_v1/report.json" \
  --native-sha256 8d75ac17b0e2d0e27a7b3151c00ac64d6f9cc959976d7a57407cc683b19fb0d6 \
  --native-representation-root "$ROOT" --output "$OUTPUT"
