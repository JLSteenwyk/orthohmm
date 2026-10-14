#!/bin/bash
#SBATCH --job-name=ohmm-boundary-native
#SBATCH --partition=gpu
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00

set -euo pipefail
source_root=$1
plan=$2
export PYTHONPATH="$source_root"
export PYTHONHASHSEED=0
export CUDA_VISIBLE_DEVICES=""
export OMP_NUM_THREADS=1
export NUMBA_NUM_THREADS=16
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
cd "$source_root"
exec /home/bizon/anaconda3/bin/python benchmark_tools/run_boundary_clear_pipeline.py run --execution-plan "$plan"
