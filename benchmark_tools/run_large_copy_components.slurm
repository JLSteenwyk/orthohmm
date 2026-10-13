#!/bin/bash
#SBATCH --job-name=ohmm-copy-components
#SBATCH --partition=gpu
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=72:00:00

set -euo pipefail
export CUDA_VISIBLE_DEVICES=""
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONHASHSEED=0 NUMBA_NUM_THREADS=16
exec /home/bizon/anaconda3/bin/python "$1/benchmark_tools/run_large_copy_components.py" run --execution-plan "$2"
