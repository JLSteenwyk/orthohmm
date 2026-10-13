#!/bin/bash
#SBATCH --job-name=ohmm-native-readback
#SBATCH --partition=gpu
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=30:00

set -euo pipefail
export CUDA_VISIBLE_DEVICES="" PYTHONFAULTHANDLER=1 PYTHONHASHSEED=0
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
exec /home/bizon/anaconda3/bin/python "$1/benchmark_tools/probe_component_native_recovery.py" --execution-plan "$2" --output "$3"
