#!/bin/bash
#SBATCH --job-name=lineage_overhead
#SBATCH --partition=spark
#SBATCH --nodelist=spark-7ff0
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=96G
#SBATCH --exclusive
#SBATCH --array=0-17%1
#SBATCH --time=01:00:00
#SBATCH --no-requeue
#SBATCH --chdir=/home/jlsteenwyk/projects/orthohmm-publication/lineage_overhead_recipe_v1
#SBATCH --output=/home/jlsteenwyk/projects/orthohmm-publication/lineage_overhead_%A_%a.log
set -euo pipefail
[[ ${1:-} =~ ^[0-9a-f]{64}$ && ${2:-} =~ ^[0-9a-f]{64}$ ]]
[[ -z ${LD_PRELOAD:-} && -z ${LD_LIBRARY_PATH:-} && -z ${LD_AUDIT:-} && -z ${PYTHONPATH:-} ]]
ROOT=/home/jlsteenwyk/projects/orthohmm-publication
export PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
export PYTHONPYCACHEPREFIX="$ROOT/lineage_collector_overhead_v1/cache_${SLURM_ARRAY_TASK_ID}"
[[ ! -e "$PYTHONPYCACHEPREFIX" ]]
exec "$ROOT/envs/orthohmm/bin/python" -B -m benchmark_tools.run_lineage_overhead_panel \
  --plan "$ROOT/lineage_overhead_recipe_v1/benchmark_tools/results/dgx_lineage_overhead_plan_20260919.json" \
  --recipe "$ROOT/lineage_overhead_recipe_v1.json" --recipe-sha "$1" \
  --authorization "$ROOT/lineage_overhead_authorization_v1.json" --authorization-sha "$2" \
  --index "${SLURM_ARRAY_TASK_ID}"
