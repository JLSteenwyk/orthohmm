#!/bin/bash
#SBATCH --job-name=ob_seq_features
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/ob_sequence_features_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONHASHSEED=0
cd "$ROOT"
exec /home/bizon/anaconda3/bin/python "$ROOT/benchmarks/work/publication_ob_sequence_features_v1/benchmark_tools/prepare_ob_sequence_features.py" --root "$ROOT" --output "$ROOT/benchmarks/results/ob_sequence_features_v1"
