#!/bin/bash
#SBATCH --job-name=ohmm-singleton-tie-ob
#SBATCH --partition=gpu
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=4-00:00:00

set -euo pipefail
source_root=$1
plan=$2
export PYTHONPATH="$source_root"
export PYTHONHASHSEED=0
export CUDA_VISIBLE_DEVICES=""
export PYTHONFAULTHANDLER=1
export OMP_NUM_THREADS=1
export NUMBA_NUM_THREADS=16
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
cd "$source_root"
exec /home/bizon/anaconda3/bin/python benchmark_tools/run_singleton_tie_orthobench.py run --execution-plan "$plan"
