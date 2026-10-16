#!/bin/bash
#SBATCH --job-name=ohmm-shared-export-recovery
#SBATCH --partition=gpu
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --no-requeue

set -euo pipefail
source_root=$1
plan=$2
export PYTHONPATH="$source_root"
export PYTHONHASHSEED=0
export CUDA_VISIBLE_DEVICES=""
export PYTHONFAULTHANDLER=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
cd "$source_root"
exec /usr/bin/time -v -o "$(dirname "$plan")/recovery.time.txt" \
    /home/bizon/anaconda3/bin/python benchmark_tools/recover_shared_qfo_export.py run --execution-plan "$plan"
