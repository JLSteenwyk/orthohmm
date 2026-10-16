#!/bin/bash
#SBATCH --job-name=ohmm-fresh-protected-union
#SBATCH --partition=gpu
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --no-requeue

set -euo pipefail
source_root=$1
plan=$2
export PYTHONPATH="$source_root"
export PYTHONHASHSEED=0
export CUDA_VISIBLE_DEVICES=""
export OMP_NUM_THREADS=8
export NUMBA_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
cd "$source_root"
exec /usr/bin/time -v -o "$(dirname "$plan")/diagnostic.time.txt" \
    /home/bizon/anaconda3/bin/python benchmark_tools/fresh_protected_union_search.py run --execution-plan "$plan"
