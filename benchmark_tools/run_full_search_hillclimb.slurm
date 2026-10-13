#!/bin/bash
#SBATCH --job-name=ohmm-full-cap
#SBATCH --partition=gpu
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --array=0-2%1

set -euo pipefail
export CUDA_VISIBLE_DEVICES=""
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONHASHSEED=0
export NUMBA_NUM_THREADS=32
datasets=(OrthoBench QfO ThreeKingdoms)
exec /home/bizon/anaconda3/bin/python "$1/benchmark_tools/run_full_search_hillclimb.py" run \
    --execution-plan "$2" --dataset "${datasets[$SLURM_ARRAY_TASK_ID]}"
