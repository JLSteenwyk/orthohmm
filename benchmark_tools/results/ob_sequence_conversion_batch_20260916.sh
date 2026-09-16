#!/bin/bash
#SBATCH --job-name=ob_seq_convert
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/ob_seq_convert_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
exec /home/bizon/anaconda3/bin/python "$ROOT/benchmarks/work/publication_ob_sequence_conversion_v1/benchmark_tools/convert_sequence_search_control.py" --manifest "$ROOT/benchmark_tools/results/ob_sequence_search_prepared_20260916.json" --output "$ROOT/benchmarks/results/ob_sequence_numeric_v1"
