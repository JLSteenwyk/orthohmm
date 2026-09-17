#!/bin/bash
#SBATCH --job-name=qfo_checked_full
#SBATCH --cpus-per-task=32
#SBATCH --mem=192G
#SBATCH --time=24:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_checked_full_replay_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
cd "$ROOT"
exec /home/bizon/anaconda3/bin/python "$ROOT/benchmarks/work/publication_qfo_checked_full_replay_v1/benchmark_tools/run_qfo_checked_full_replay.py" --root "$ROOT" --output "$ROOT/benchmarks/results/qfo_checked_full_replay_v1"
