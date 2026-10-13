#!/bin/bash
#SBATCH --job-name=ohmm-component-admit
#SBATCH --partition=gpu
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=01:00:00

set -euo pipefail
export CUDA_VISIBLE_DEVICES="" PYTHONFAULTHANDLER=1 PYTHONHASHSEED=0
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=16
exec /home/bizon/anaconda3/bin/python "$1/benchmark_tools/admit_component_recovery.py" --original-plan "$2" --candidate-report "$3" --scoring-report "$4" --output "$5"
