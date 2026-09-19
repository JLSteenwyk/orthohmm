#!/bin/bash
#SBATCH --job-name=dual_native
#SBATCH --partition=spark
#SBATCH --nodelist=spark-7ff0
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=96G
#SBATCH --exclusive
#SBATCH --time=01:00:00
#SBATCH --no-requeue
#SBATCH --output=/home/jlsteenwyk/projects/orthohmm-publication/dual_native_%j.log
set -euo pipefail
[[ ${1:-} =~ ^[0-9a-f]{64}$ ]]
[[ ${2:-} =~ ^[012]$ ]]
cd /home/jlsteenwyk/projects/orthohmm-publication/dual_native_recipe_v1
export PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1
exec /home/jlsteenwyk/projects/orthohmm-publication/envs/orthohmm/bin/python -B \
  -m benchmark_tools.run_dual_native_diagnostic \
  --plan benchmark_tools/results/dgx_dual_native_plan_20260919.json \
  --plan-sha 029f19b0e21f356387ae4e2df50310a225cdba4f356037814646af57226e1113 \
  --recipe /home/jlsteenwyk/projects/orthohmm-publication/dual_native_recipe_v1.json \
  --recipe-sha "$1" --index "$2"
