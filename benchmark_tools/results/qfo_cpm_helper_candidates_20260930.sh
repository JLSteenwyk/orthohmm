#!/bin/bash
#SBATCH --job-name=qfo_cpm_helper_candidates
#SBATCH --nodelist=bizon
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_cpm_helper_candidates_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Frozen executor required}
COMMIT=${2:?Exact executor commit required}
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --quiet HEAD -- benchmark_tools orthohmm
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
unset PYTHONPATH PYTHONHOME LD_PRELOAD LD_LIBRARY_PATH
cd "$ROOT"
exec /usr/bin/time -v -o "$ROOT/benchmarks/work/qfo_cpm_helper_candidates_${SLURM_JOB_ID}.time.txt" \
  /home/bizon/anaconda3/bin/python -B "$EXECUTOR/benchmark_tools/prepare_recovered_cpm_candidates.py" \
  --root "$ROOT" \
  --helper-recovery-readback "$ROOT/benchmark_tools/results/qfo_cpm_helper_recovery_admission_readback_20260930.json" \
  --helper-recovery-readback-sha256 5fb31ea6cee69501d5ad1826789e9a178f230bb909012537e123d18d21d65850 \
  --helper-candidate-protocol-sha256 f07242533466ac094bb2151117b1289169a05b6ecface05e4ff9655216ef8780
