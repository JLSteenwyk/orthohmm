#!/bin/bash
#SBATCH --job-name=qfo_replay_cap
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_replay_capture_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONHASHSEED=0
exec /home/bizon/anaconda3/bin/python "$ROOT/benchmarks/work/publication_qfo_replay_capture_v1/benchmark_tools/capture_qfo_replay_graph.py" --root "$ROOT" --output "$ROOT/benchmarks/results/qfo_replay_initial_capture_v1"
