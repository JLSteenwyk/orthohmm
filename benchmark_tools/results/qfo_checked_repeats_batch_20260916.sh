#!/bin/bash
#SBATCH --job-name=qfo_checked_repeat
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_checked_repeats_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
cd "$ROOT"
exec /home/bizon/anaconda3/bin/python "$ROOT/benchmarks/work/publication_qfo_checked_repeats_v1/benchmark_tools/run_qfo_checked_repeats.py" --root "$ROOT" --output "$ROOT/benchmarks/results/qfo_checked_repeats_v1"
