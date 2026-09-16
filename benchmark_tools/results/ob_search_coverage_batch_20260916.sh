#!/bin/bash
#SBATCH --job-name=ob_hit_coverage
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/ob_hit_coverage_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONHASHSEED=0
exec /home/bizon/anaconda3/bin/python "$ROOT/benchmarks/work/publication_ob_hit_coverage_v1/benchmark_tools/compare_search_hit_coverage.py" --root "$ROOT" --output "$ROOT/benchmarks/results/ob_search_hit_coverage_v1.json"
