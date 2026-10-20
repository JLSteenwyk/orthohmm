#!/bin/bash
#SBATCH --job-name=ohmm-weighted-repair-selection
#SBATCH --partition=gpu
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00

set -euo pipefail
export CUDA_VISIBLE_DEVICES=""
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export NUMBA_NUM_THREADS=1 PYTHONHASHSEED=0
exec /home/bizon/anaconda3/bin/python "$1/benchmark_tools/plan_weighted_projection_followup.py" \
    --repo "$2" --driver "$1" --output "$3"
