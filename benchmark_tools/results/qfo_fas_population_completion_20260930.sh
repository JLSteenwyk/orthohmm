#!/bin/bash
#SBATCH --job-name=qfo-fas-completion
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --no-requeue

set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
/home/bizon/anaconda3/bin/python -B -m benchmark_tools.complete_fas_population \
  --output benchmarks/work/qfo_fas_population_completion_20260930
