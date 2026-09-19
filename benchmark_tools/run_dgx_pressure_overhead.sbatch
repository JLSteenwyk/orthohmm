#!/bin/bash
#SBATCH --job-name=dgx_pressure_overhead
#SBATCH --partition=spark
#SBATCH --nodelist=spark-7ff0
#SBATCH --exclusive
#SBATCH --cpus-per-task=20
#SBATCH --mem=96G
#SBATCH --array=0-17%1
#SBATCH --time=01:00:00
#SBATCH --no-requeue
#SBATCH --chdir=/home/jlsteenwyk/projects/orthohmm-publication/pressure_overhead_recipe_v1
#SBATCH --output=/home/jlsteenwyk/projects/orthohmm-publication/pressure_overhead_%A_%a.log

set -euo pipefail
ROOT=/home/jlsteenwyk/projects/orthohmm-publication
exec env -u LD_PRELOAD -u LD_LIBRARY_PATH -u LD_AUDIT PYTHONHASHSEED=0 \
  /usr/bin/python3 -B "$ROOT/pressure_overhead_recipe_v1/benchmark_tools/run_dgx_frontier_overhead.py" \
  --plan "$ROOT/pressure_overhead_recipe_v1/benchmark_tools/results/dgx_pressure_overhead_plan_20260919.json" \
  --plan-sha 3950ccfa867c463ccfd7d8693dc331a3c85dc75ba67d9e8db89be4c1dda91c38 \
  --recipe "$ROOT/pressure_overhead_recipe_v1_manifest.json" \
  --recipe-sha 11ce57a7895d92c456cc2b438157355cb0773cd473f4fbc4fc8305330a33eda3 \
  --authorization "$ROOT/pressure_overhead_authorization_v1.json" \
  --authorization-sha d1f5878580ed85829b1ce50ba0eebbcfd573975a6eb2b28471762a81978b9dea \
  --index "${SLURM_ARRAY_TASK_ID}"
