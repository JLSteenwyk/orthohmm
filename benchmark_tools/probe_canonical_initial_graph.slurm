#!/bin/bash
#SBATCH --job-name=ohmm-canonical-graph-probe
#SBATCH --partition=gpu
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=0-06:00:00

set -euo pipefail
source_root=$1
plan=$2
export PYTHONHASHSEED=0
export CUDA_VISIBLE_DEVICES=""
export PYTHONFAULTHANDLER=1
export OMP_NUM_THREADS=1
export NUMBA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
cd "$source_root"
exec /home/bizon/anaconda3/bin/python benchmark_tools/probe_canonical_initial_graph.py run --execution-plan "$plan"
