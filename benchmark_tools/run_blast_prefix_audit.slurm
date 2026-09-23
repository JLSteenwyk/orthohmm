#!/bin/bash
#SBATCH --job-name=qfo_blast_prefix_audit
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_blast_prefix_audit_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Pinned executor required}
COMMIT=${2:?Exact revision required}
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --exit-code HEAD -- benchmark_tools
cd "$ROOT"
unset LD_PRELOAD LD_LIBRARY_PATH LD_AUDIT PYTHONPATH
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
exec /home/bizon/anaconda3/bin/python "$EXECUTOR/benchmark_tools/audit_interrupted_blast_rows.py" \
    --root "$ROOT" --output "$ROOT/benchmarks/work/qfo_blast_prefix_audit_20260923"
