#!/bin/bash
#SBATCH --job-name=qfo_recovered_reference_impact
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_recovered_reference_impact_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR="$ROOT/benchmarks/work/publication_native_failure_membership_v1_20260925"
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == 9b9ab8dc7f56b7dfc29e2eec2f94cc919b06a018 ]]
[[ -z $(git -C "$EXECUTOR" status --porcelain) ]]
AUDIT="$ROOT/benchmark_tools/results/qfo_recovered_failure_membership_20260926.json"
[[ $(sha256sum "$AUDIT" | cut -d ' ' -f 1) == bf7093d82dff31158fe5474d9576603fef3fcc83a62b7761e27eaa1355cbc79e ]]
OUTPUT="$ROOT/benchmarks/work/qfo_recovered_reference_impact_20260926"
[[ ! -e "$OUTPUT" && ! -L "$OUTPUT" ]]
cd "$ROOT"
exec env -i HOME="$HOME" USER=bizon LOGNAME=bizon PATH=/usr/local/bin:/usr/bin:/bin LANG=C LC_ALL=C \
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  "$ROOT/benchmarks/work/orthomcl_python_env_20260918/bin/python" -I -B \
  "$EXECUTOR/benchmark_tools/audit_orthomcl_reference_impact.py" \
  --audit "$AUDIT" \
  --reference "$ROOT/qfo_benchmark/benchmark-webservice/reference_data/2020" \
  --benchmark-repo "$ROOT/qfo_benchmark/benchmark-webservice" \
  --darwin-image "$ROOT/qfo_benchmark/scoring/container_cache/qfobenchmark-darwin-2022.1.img" \
  --output "$OUTPUT"
