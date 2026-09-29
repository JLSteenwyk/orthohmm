#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --no-requeue
#SBATCH --dependency=afterany:22377:22378
#SBATCH --job-name=oh-completion-audits
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/publication_completion_audits_20260929_%j.log
set -euo pipefail
cd /mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
exec env -u LD_PRELOAD -u LD_LIBRARY_PATH -u LD_AUDIT \
  /home/bizon/anaconda3/bin/python -B -m benchmark_tools.run_publication_completion_audits "$@"
