#!/bin/bash
#SBATCH --job-name=qfo-fas-population
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:30:00
#SBATCH --no-requeue

set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
/home/bizon/anaconda3/bin/python -B -m benchmark_tools.audit_fas_population \
  --output benchmarks/work/qfo_fas_population_20260930
