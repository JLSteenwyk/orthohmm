#!/bin/bash
#SBATCH --job-name=ob_tree_score
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/ob_tree_scoring_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONHASHSEED=0
cd "$ROOT"
exec /home/bizon/anaconda3/bin/python "$ROOT/benchmarks/work/publication_ob_tree_scoring_v1/benchmark_tools/assemble_species_tree_robustness.py" --root "$ROOT" --output "$ROOT/benchmarks/results/ob_species_tree_robustness_scoring_v1"
