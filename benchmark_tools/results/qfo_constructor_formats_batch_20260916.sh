#!/bin/bash
#SBATCH --job-name=qfo_edge_formats
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_constructor_formats_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
cd "$ROOT"
exec /home/bizon/anaconda3/bin/python "$ROOT/benchmarks/work/publication_qfo_constructor_formats_v1/benchmark_tools/probe_qfo_direct_graph.py" --root "$ROOT" --output "$ROOT/benchmarks/results/qfo_constructor_formats_v1" --compare-formats
