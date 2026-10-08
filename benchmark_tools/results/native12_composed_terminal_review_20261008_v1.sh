#!/bin/bash
#SBATCH --job-name=ohmm_native12_review
#SBATCH --partition=gpu
#SBATCH --nodelist=bizon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=128G
#SBATCH --time=06:00:00
#SBATCH --no-requeue
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
[[ "${SLURM_JOB_ID:?}" =~ ^[0-9]+$ && "${SLURM_CPUS_PER_TASK:?}" == 2 && "${SLURM_MEM_PER_NODE:?}" == 131072 ]]
[[ "$#" == 3 ]]
cd "$ROOT"
unset PYTHONPATH PYTHONHOME PYTHONUSERBASE LD_PRELOAD LD_LIBRARY_PATH LD_AUDIT
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
exec "$ROOT/benchmarks/work/native_factorial_review_py310_20261004/bin/python" -B -X faulthandler \
    -m benchmark_tools.review_native12_composed_attempt --request "$1" --request-sha256 "$2" --source-sha256 "$3"
