#!/bin/bash
#SBATCH --job-name=root_context_overhead
#SBATCH --partition=spark
#SBATCH --nodelist=spark-7ff0
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=96G
#SBATCH --exclusive
#SBATCH --time=05:00:00
#SBATCH --no-requeue
#SBATCH --output=/home/jlsteenwyk/projects/orthohmm-publication/root_context_overhead_%j.log
set -euo pipefail
[[ ${1:-} =~ ^[0-9a-f]{64}$ ]]
[[ -z ${LD_PRELOAD:-} && -z ${LD_LIBRARY_PATH:-} && -z ${LD_AUDIT:-} && -z ${PYTHONPATH:-} ]]
cd /home/jlsteenwyk/projects/orthohmm-publication/root_context_overhead_recipe_v1
export PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
export PYTHONPYCACHEPREFIX=/home/jlsteenwyk/projects/orthohmm-publication/root_context_overhead_v1/cache_observer
[[ ! -e "$PYTHONPYCACHEPREFIX" ]]
exec /home/jlsteenwyk/projects/orthohmm-publication/envs/orthohmm/bin/python -B \
  -m benchmark_tools.run_root_context_overhead \
  --plan benchmark_tools/results/dgx_root_context_overhead_plan_20260919.json \
  --recipe /home/jlsteenwyk/projects/orthohmm-publication/root_context_overhead_recipe_v1.json \
  --recipe-sha "$1"
