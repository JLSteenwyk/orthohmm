#!/bin/bash
#SBATCH --job-name=swiss_identity_inputs
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/corrected_swiss_identity_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR="$ROOT/benchmarks/work/publication_swiss_identity_v1"
test "$(git -C "$EXECUTOR" rev-parse HEAD)" = "$(git -C "$ROOT" rev-parse '21cf7bb^{commit}')"
git -C "$EXECUTOR" diff --exit-code HEAD -- benchmark_tools
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONHASHSEED=0
cd "$EXECUTOR"
exec /home/bizon/anaconda3/bin/python -m benchmark_tools.prepare_corrected_swiss_alignments \
  --root "$ROOT" --output "$ROOT/benchmarks/results/corrected_swiss_identity_v1" \
  --protocol-sha256 aaaeabebb2cb3ec5e95dc38f99b53f6814579e1bbdf7241da5e4391be6a35f55
