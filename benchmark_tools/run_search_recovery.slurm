#!/bin/bash
#SBATCH --job-name=ohmm-hmm-recovery
#SBATCH --partition=gpu
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --array=0-5%1

set -euo pipefail
source_root=$1
plan=$2
arms=(control cap400 shared_hit1 reduced_k4 exact_k3 shortlist20000)
export PYTHONPATH="$source_root"
export PYTHONHASHSEED=0
export CUDA_VISIBLE_DEVICES=""
export OMP_NUM_THREADS=8
export NUMBA_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
cd "$source_root"
/home/bizon/anaconda3/bin/python benchmark_tools/run_search_recovery.py run \
    --execution-plan "$plan" --label "${arms[$SLURM_ARRAY_TASK_ID]}"
/home/bizon/anaconda3/bin/python benchmark_tools/audit_search_recovery.py \
    --execution-plan "$plan" --output "$(dirname "$plan")/summary.json"
