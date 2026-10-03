#!/bin/bash
#SBATCH --job-name=orthohmm_shared_scaling
#SBATCH --nodelist=bizon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --time=1-02:00:00
#SBATCH --no-requeue
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
unset PYTHONHOME LD_PRELOAD LD_LIBRARY_PATH LD_AUDIT
PYTHON="$ROOT/benchmarks/work/threadripper_private_controller_20260928/venv/bin/python"
read -r ACTUAL _ < <(sha256sum -- "$PYTHON")
[[ "$ACTUAL" == 8b1cd756be711ef53f35cb6c954472fdfc52094c4619a01f553f64354587388b ]]
cd "$ROOT"
export ORTHOHMM_THREADRIPPER_DEPLOYMENT=private_v2_20260928
export PYTHONHASHSEED=0 PYTHONNOUSERSITE=1
export PYTHONPYCACHEPREFIX=/dev/shm/orthohmm_scaling_driver_${SLURM_JOB_ID:?}
[[ ! -e "$PYTHONPYCACHEPREFIX" && ! -L "$PYTHONPYCACHEPREFIX" ]]
exec "$PYTHON" -B benchmark_tools/run_threadripper_scaling.py "$@"
