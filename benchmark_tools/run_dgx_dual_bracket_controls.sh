#!/bin/bash
#SBATCH --job-name=dual_cpu_controls
#SBATCH --partition=spark
#SBATCH --nodelist=spark-7ff0
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=256M
#SBATCH --exclusive
#SBATCH --time=00:05:00
#SBATCH --no-requeue
#SBATCH --output=/home/jlsteenwyk/projects/orthohmm-publication/dual_bracket_controls_%j.log
set -euo pipefail
cd /home/jlsteenwyk/projects/orthohmm-publication/dual_bracket_recipe_v1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
exec /home/jlsteenwyk/projects/orthohmm-publication/envs/orthohmm/bin/python -B \
  -m benchmark_tools.run_dual_bracket_controls \
  --output /home/jlsteenwyk/projects/orthohmm-publication/dual_bracket_controls_v1
