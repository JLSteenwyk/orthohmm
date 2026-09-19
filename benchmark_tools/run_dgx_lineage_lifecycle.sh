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
#SBATCH --output=/home/jlsteenwyk/projects/orthohmm-publication/lineage_lifecycle_%j.log
set -euo pipefail
cd /home/jlsteenwyk/projects/orthohmm-publication/lineage_lifecycle_source_v1
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
exec /usr/bin/python3 -B -m benchmark_tools.run_lineage_lifecycle_control \
  --output /home/jlsteenwyk/projects/orthohmm-publication/lineage_lifecycle_v1
