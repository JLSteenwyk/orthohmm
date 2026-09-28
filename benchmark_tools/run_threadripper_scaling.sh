#!/bin/bash
#SBATCH --job-name=orthohmm_local_scaling
#SBATCH --nodelist=bizon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --exclusive
#SBATCH --time=1-02:00:00
#SBATCH --no-requeue
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
cd "$ROOT"
export PYTHONHASHSEED=0 PYTHONNOUSERSITE=1
export PYTHONPYCACHEPREFIX=/dev/shm/orthohmm_scaling_driver_${SLURM_JOB_ID:?}
[[ ! -e "$PYTHONPYCACHEPREFIX" && ! -L "$PYTHONPYCACHEPREFIX" ]]
unset PYTHONHOME LD_PRELOAD LD_LIBRARY_PATH LD_AUDIT
exec /home/bizon/anaconda3/bin/python -B benchmark_tools/run_threadripper_scaling.py "$@"
