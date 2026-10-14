#!/bin/bash
#SBATCH --job-name=ohmm-singleton-followup
#SBATCH --partition=gpu
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=10-00:00:00

set -euo pipefail
driver=$1
plan=$2
export PYTHONHASHSEED=0
export CUDA_VISIBLE_DEVICES=""
export PYTHONFAULTHANDLER=1
export OMP_NUM_THREADS=1
export NUMBA_NUM_THREADS=16
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
exec /home/bizon/anaconda3/bin/python "$driver/benchmark_tools/run_singleton_tie_followup.py" run --execution-plan "$plan"
