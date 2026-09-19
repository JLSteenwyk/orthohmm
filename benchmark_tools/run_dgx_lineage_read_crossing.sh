#!/bin/bash
#SBATCH --job-name=lineage_read_crossing
#SBATCH --partition=spark
#SBATCH --nodelist=spark-7ff0
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --exclusive
#SBATCH --time=00:10:00
#SBATCH --no-requeue
#SBATCH --chdir=/home/jlsteenwyk/projects/orthohmm-publication
#SBATCH --export=ALL,TMPDIR=/tmp
#SBATCH --output=/home/jlsteenwyk/projects/orthohmm-publication/lineage_read_crossing_%j.log
set -euo pipefail
cd /home/jlsteenwyk/projects/orthohmm-publication/lineage_read_crossing_source_v1
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
exec srun --exclusive --exact --nodes=1 --ntasks=1 --cpus-per-task=2 --cpu-bind=cores \
  /usr/bin/python3 -B -m benchmark_tools.run_lineage_read_crossing_control \
  --output /home/jlsteenwyk/projects/orthohmm-publication/lineage_read_crossing_v1
