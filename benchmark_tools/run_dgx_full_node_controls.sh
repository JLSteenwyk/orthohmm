#!/bin/bash
#SBATCH --job-name=full_node_controls
#SBATCH --partition=spark
#SBATCH --nodelist=spark-7ff0
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=96G
#SBATCH --exclusive
#SBATCH --time=00:45:00
#SBATCH --no-requeue
#SBATCH --output=/home/jlsteenwyk/projects/orthohmm-publication/full_node_controls_%j.log
set -euo pipefail
[[ ${1:-} =~ ^[0-9a-f]{64}$ ]]
[[ -z ${LD_PRELOAD:-} && -z ${LD_LIBRARY_PATH:-} && -z ${LD_AUDIT:-} && -z ${PYTHONPATH:-} ]]
cd /home/jlsteenwyk/projects/orthohmm-publication/full_node_controls_recipe_v1
export PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
export PYTHONPYCACHEPREFIX=/home/jlsteenwyk/projects/orthohmm-publication/full_node_controls_cache_v1
[[ ! -e "$PYTHONPYCACHEPREFIX" ]]
exec /home/jlsteenwyk/projects/orthohmm-publication/envs/orthohmm/bin/python -B \
  -m benchmark_tools.run_full_node_controls \
  --output /home/jlsteenwyk/projects/orthohmm-publication/full_node_controls_v1 \
  --recipe /home/jlsteenwyk/projects/orthohmm-publication/full_node_controls_recipe_v1.json \
  --recipe-sha "$1"
