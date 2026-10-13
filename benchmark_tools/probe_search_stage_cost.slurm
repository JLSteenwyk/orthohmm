#!/bin/bash
#SBATCH --job-name=ohmm-search-stages
#SBATCH --partition=gpu
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=20:00

set -euo pipefail
export CUDA_VISIBLE_DEVICES=""
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export NUMBA_NUM_THREADS=4 PYTHONHASHSEED=0 PYTHONFAULTHANDLER=1
exec /home/bizon/anaconda3/bin/python "$1/benchmark_tools/probe_search_stage_cost.py" --repo "$2" --source-root "$1" --output "$3"
