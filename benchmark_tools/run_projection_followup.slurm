#!/bin/bash
#SBATCH --job-name=ohmm-projection-transfer
#SBATCH --partition=gpu
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=1-00:00:00

set -euo pipefail
export CUDA_VISIBLE_DEVICES=""
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export NUMBA_NUM_THREADS=16 PYTHONHASHSEED=0
exec /home/bizon/anaconda3/bin/python "$1/benchmark_tools/run_projection_followup.py" run \
    --execution-plan "$2"
