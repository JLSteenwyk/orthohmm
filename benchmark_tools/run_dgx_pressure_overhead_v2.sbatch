#!/bin/bash
#SBATCH --job-name=dgx_pressure_overhead_v2
#SBATCH --partition=spark
#SBATCH --nodelist=spark-7ff0
#SBATCH --exclusive
#SBATCH --cpus-per-task=20
#SBATCH --mem=96G
#SBATCH --array=0-17%1
#SBATCH --time=01:00:00
#SBATCH --no-requeue
#SBATCH --chdir=/home/jlsteenwyk/projects/orthohmm-publication/pressure_overhead_recipe_v2
#SBATCH --output=/home/jlsteenwyk/projects/orthohmm-publication/pressure_overhead_v2_%A_%a.log

set -euo pipefail
ROOT=/home/jlsteenwyk/projects/orthohmm-publication
exec env -u LD_PRELOAD -u LD_LIBRARY_PATH -u LD_AUDIT PYTHONHASHSEED=0 \
  "$ROOT/envs/orthohmm/bin/python" -B "$ROOT/pressure_overhead_recipe_v2/benchmark_tools/run_dgx_frontier_overhead.py" \
  --plan "$ROOT/pressure_overhead_recipe_v2/benchmark_tools/results/dgx_pressure_overhead_plan_v2_20260919.json" \
  --plan-sha b644e165dbf4d0beabf1cf4d9b6c314de522e3ebd1b91598ebebea99094c8fff \
  --recipe "$ROOT/pressure_overhead_recipe_v2_manifest.json" \
  --recipe-sha 50fc3a14c4d5b53c4eb158a22efcfb8b289b346718fad3c78b300172e319d8a1 \
  --authorization "$ROOT/pressure_overhead_authorization_v2.json" \
  --authorization-sha 17943df08dfe21bf1975ff0dcb29dc66a9507b156e4698b1672e7517cdfd85e6 \
  --index "${SLURM_ARRAY_TASK_ID}"
