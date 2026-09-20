#!/bin/bash
#SBATCH --job-name=root_context_scaling
#SBATCH --partition=spark
#SBATCH --nodelist=spark-7ff0
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=96G
#SBATCH --exclusive
#SBATCH --time=1-00:00:00
#SBATCH --no-requeue
#SBATCH --output=/home/jlsteenwyk/projects/orthohmm-publication/root_context_scaling_%j.log
set -euo pipefail
[[ $# -eq 4 ]]
[[ $1 =~ ^([0-9]|1[0-9]|2[0-6])$ ]]
[[ $2 =~ ^[0-9a-f]{64}$ && $4 =~ ^[0-9a-f]{64}$ ]]
[[ $3 = /* ]]
[[ -z ${LD_PRELOAD:-} && -z ${LD_LIBRARY_PATH:-} && -z ${LD_AUDIT:-} && -z ${PYTHONPATH:-} ]]
cd /home/jlsteenwyk/projects/orthohmm-publication/root_context_scaling_recipe_v1
export PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
printf -v suffix '%02d' "$1"
export PYTHONPYCACHEPREFIX="/home/jlsteenwyk/projects/orthohmm-publication/scaling_root_context_v1/cache_observer_${suffix}"
[[ ! -e "$PYTHONPYCACHEPREFIX" && ! -L "$PYTHONPYCACHEPREFIX" ]]
exec /home/jlsteenwyk/projects/orthohmm-publication/envs/orthohmm/bin/python -B \
  -m benchmark_tools.run_root_context_scaling \
  --plan benchmark_tools/results/dgx_root_context_scaling_plan_v2_20260920.json \
  --recipe /home/jlsteenwyk/projects/orthohmm-publication/root_context_scaling_recipe_v1.json \
  --index "$1" --recipe-sha "$2" --authorization "$3" --authorization-sha "$4"
