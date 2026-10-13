#!/bin/bash
#SBATCH --job-name=ohmm-membership-r1
#SBATCH --partition=gpu
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=2-00:00:00
#SBATCH --array=0-2%1

set -euo pipefail

source_root=$1
plan=$2
datasets=(OrthoBench QfO ThreeKingdoms)
export PYTHONPATH="$source_root"
export PYTHONHASHSEED=0
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

cd "$source_root"
exec /home/bizon/anaconda3/bin/python benchmark_tools/run_membership_hillclimb.py run \
    --plan "$plan" --dataset "${datasets[$SLURM_ARRAY_TASK_ID]}"
