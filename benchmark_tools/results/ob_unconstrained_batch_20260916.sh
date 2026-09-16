#!/bin/bash
#SBATCH --job-name=ob_unconstrained
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/ob_unconstrained_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
cd "$ROOT"
exec /home/bizon/anaconda3/bin/python "$ROOT/benchmarks/work/publication_ob_unconstrained_executor_v1/benchmark_tools/run_orthobench_factorial_cell.py" \
  --manifest "$ROOT/benchmark_tools/results/orthobench_factorial_prepared_20260916.json" \
  --manifest-sha256 5c325f4d77865e0c7571fe4bb4df0be0959977a49e1d22169f192518700f9382 \
  --environment-manifest "$ROOT/benchmark_tools/results/publication_variable_native_methods_20260916.json" \
  --environment-sha256 bf677728cf9312edd9755d9ac1ceefbce3abb8c3848d43413631dc00ca00eb2f \
  --preparation-job 21161 --index 3 --unconstrained
