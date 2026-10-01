#!/bin/bash
#SBATCH --job-name=qfo_private_phylogeny_admission
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_private_phylogeny_admission_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Frozen admission executor required}
COMMIT=${2:?Exact executor commit required}
PROTOCOL_SHA=${3:?Reviewed admission protocol required}
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --quiet HEAD -- benchmark_tools orthohmm
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
unset PYTHONNOUSERSITE PYTHONPATH PYTHONHOME LD_PRELOAD LD_LIBRARY_PATH LD_AUDIT
cd "$ROOT"
exec /home/bizon/anaconda3/bin/python -B "$EXECUTOR/benchmark_tools/admit_qfo_private_phylogeny_control.py" \
  --root "$ROOT" --protocol-sha256 "$PROTOCOL_SHA" \
  --submission-sha256 5eeafb0d9c478c06c5e110fb0b6a8638bd1d958e3ad278808b824cfb4028d9e3 \
  --output "$ROOT/benchmarks/work/qfo_private_phylogeny_control_admission_20261001.json"
