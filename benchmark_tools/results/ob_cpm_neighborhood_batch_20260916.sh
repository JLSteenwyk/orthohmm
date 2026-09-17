#!/bin/bash
#SBATCH --job-name=ob_cpm_neighborhood
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/ob_cpm_neighborhood_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
cd "$ROOT"
exec /home/bizon/anaconda3/bin/python "$ROOT/benchmarks/work/publication_ob_cpm_neighborhood_v1/benchmark_tools/run_ob_cpm_neighborhood.py" --root "$ROOT" --output "$ROOT/benchmarks/results/ob_cpm_neighborhood_v1"
