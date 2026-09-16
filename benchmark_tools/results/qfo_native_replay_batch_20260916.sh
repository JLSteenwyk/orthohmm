#!/bin/bash
#SBATCH --job-name=qfo_native_replay
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_native_replay_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
exec /home/bizon/anaconda3/bin/python "$ROOT/benchmarks/work/publication_qfo_replay_executor_v1/benchmark_tools/run_qfo_publication_replay.py" --root "$ROOT" --output "$ROOT/benchmarks/results/publication_qfo_replay_check_v1"
