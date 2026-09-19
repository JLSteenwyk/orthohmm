#!/bin/bash
#SBATCH --job-name=lineage_native
#SBATCH --partition=spark
#SBATCH --nodelist=spark-7ff0
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=96G
#SBATCH --exclusive
#SBATCH --time=01:00:00
#SBATCH --no-requeue
#SBATCH --chdir=/home/jlsteenwyk/projects/orthohmm-publication
#SBATCH --export=ALL,TMPDIR=/tmp
#SBATCH --output=/home/jlsteenwyk/projects/orthohmm-publication/lineage_native_%j.log
set -euo pipefail
[[ ${1:-} =~ ^[0-9a-f]{64}$ ]]
[[ ${2:-} =~ ^[012]$ ]]
cd /home/jlsteenwyk/projects/orthohmm-publication/lineage_native_recipe_v1
export PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1
export PYTHONPYCACHEPREFIX=/home/jlsteenwyk/projects/orthohmm-publication/lineage_native_v1/cache_$2
[[ ! -e "$PYTHONPYCACHEPREFIX" ]]
exec /home/jlsteenwyk/projects/orthohmm-publication/envs/orthohmm/bin/python -B \
  -m benchmark_tools.run_lineage_native_diagnostic \
  --plan benchmark_tools/results/dgx_lineage_native_plan_20260919.json \
  --plan-sha 714ef04458904e1c01526d171d0b2fde658ac135bd416af49ab107dc5ed14bfe \
  --recipe /home/jlsteenwyk/projects/orthohmm-publication/lineage_native_recipe_v1.json \
  --recipe-sha "$1" --index "$2"
