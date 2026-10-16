#!/bin/bash
#SBATCH --job-name=ohmm-protected-union
#SBATCH --partition=gpu
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --no-requeue

set -euo pipefail
source_root=$1
repo=$2
output=$3
export PYTHONPATH="$source_root"
export PYTHONHASHSEED=0
export CUDA_VISIBLE_DEVICES=""
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
cd "$source_root"
exec /usr/bin/time -v -o "${output}.time.txt" \
    /home/bizon/anaconda3/bin/python benchmark_tools/diagnose_protected_candidate_union.py \
    --repo "$repo" --output "$output"
