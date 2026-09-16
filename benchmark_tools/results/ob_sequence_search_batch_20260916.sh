#!/bin/bash
#SBATCH --job-name=ob_sequence_search
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/ob_sequence_search_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
exec /home/bizon/anaconda3/bin/python "$ROOT/benchmarks/work/publication_ob_sequence_executor_v1/benchmark_tools/run_sequence_search_control.py" --manifest "$ROOT/benchmark_tools/results/ob_sequence_search_prepared_20260916.json"
