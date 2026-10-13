#!/bin/bash
#SBATCH --job-name=ohmm-matched-k4-costs
#SBATCH --partition=gpu
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=20:00

set -euo pipefail
export CUDA_VISIBLE_DEVICES="" PYTHONFAULTHANDLER=1 PYTHONHASHSEED=0
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=4
exec /home/bizon/anaconda3/bin/python "$1/benchmark_tools/probe_matched_live_k4_search.py" --repo "$2" --output "$3"
