#!/bin/bash
#SBATCH --job-name=ohmm-deferral-trace
#SBATCH --partition=gpu
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00

set -euo pipefail
export CUDA_VISIBLE_DEVICES=""
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
exec /home/bizon/anaconda3/bin/python "$1/benchmark_tools/trace_deferred_reference_pairs.py" run --execution-plan "$2"
