#!/bin/bash
#SBATCH --job-name=ohmm-gene-lca-diag
#SBATCH --partition=gpu
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00

set -euo pipefail
export CUDA_VISIBLE_DEVICES="" PYTHONHASHSEED=0 PYTHONFAULTHANDLER=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
exec /home/bizon/anaconda3/bin/python "$1/benchmark_tools/diagnose_component_gene_lca.py" --repo "$2" --output "$3" --source-revision "$4"
