#!/bin/bash
#SBATCH --job-name=ohmm-candidate-readback
#SBATCH --partition=gpu
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=01:00:00

set -euo pipefail
export CUDA_VISIBLE_DEVICES="" PYTHONFAULTHANDLER=1 PYTHONHASHSEED=0
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=16
exec /home/bizon/anaconda3/bin/python "$1/benchmark_tools/reconstruct_component_recovery.py" --original-plan "$2" --output "$3"
