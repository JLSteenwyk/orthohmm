#!/bin/bash
#SBATCH --job-name=ohmm-component-scores
#SBATCH --partition=gpu
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00

set -euo pipefail
export CUDA_VISIBLE_DEVICES="" PYTHONFAULTHANDLER=1 PYTHONHASHSEED=0
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=16
exec /usr/bin/time -v -o "$2/recovery.time.txt" /home/bizon/anaconda3/bin/python "$1/benchmark_tools/recover_component_scores.py" run --execution-plan "$2/execution_plan.json"
