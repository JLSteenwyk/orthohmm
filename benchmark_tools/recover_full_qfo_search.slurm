#!/bin/bash
#SBATCH --job-name=ohmm-qfo-cap-recover
#SBATCH --partition=gpu
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=240:00:00

set -euo pipefail
export CUDA_VISIBLE_DEVICES=""
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export NUMBA_NUM_THREADS=32 PYTHONHASHSEED=0
exec /home/bizon/anaconda3/bin/python "$1/benchmark_tools/recover_full_qfo_search.py" run --execution-plan "$2"
