#!/bin/bash
#SBATCH --job-name=lineage_lifecycle
#SBATCH --partition=spark
#SBATCH --nodelist=spark-7ff0
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --exclusive
#SBATCH --time=00:05:00
#SBATCH --no-requeue
#SBATCH --chdir=/home/jlsteenwyk/projects/orthohmm-publication
#SBATCH --export=ALL,TMPDIR=/tmp
#SBATCH --output=/home/jlsteenwyk/projects/orthohmm-publication/lineage_lifecycle_%j.log
set -euo pipefail
cd /home/jlsteenwyk/projects/orthohmm-publication/lineage_lifecycle_source_v2
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
exec srun --exclusive --exact --nodes=1 --ntasks=1 --cpus-per-task=2 --cpu-bind=cores \
  /usr/bin/python3 -B -m benchmark_tools.run_lineage_lifecycle_control \
  --output /home/jlsteenwyk/projects/orthohmm-publication/lineage_lifecycle_v2
