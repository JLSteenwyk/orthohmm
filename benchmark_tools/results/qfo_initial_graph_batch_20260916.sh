#!/bin/bash
#SBATCH --job-name=qfo_graph_diag
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_graph_diag_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONHASHSEED=0
exec /home/bizon/anaconda3/bin/python "$ROOT/benchmarks/work/publication_qfo_graph_diagnostic_v1/benchmark_tools/diagnose_qfo_initial_graph.py" --root "$ROOT" --output "$ROOT/benchmarks/results/qfo_initial_graph_diagnostic_v1"
