#!/bin/bash
#SBATCH --job-name=orthohmm_allocated_factorial
#SBATCH --partition=gpu
#SBATCH --nodelist=bizon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --time=1-02:00:00
#SBATCH --no-requeue
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
unset PYTHONHOME PYTHONPATH PYTHONUSERBASE LD_PRELOAD LD_LIBRARY_PATH LD_AUDIT
PYTHON="$ROOT/benchmarks/work/threadripper_private_controller_20260928/venv/bin/python"
read -r ACTUAL _ < <(sha256sum -- "$PYTHON")
[[ "$ACTUAL" == 8b1cd756be711ef53f35cb6c954472fdfc52094c4619a01f553f64354587388b ]]
cd "$ROOT"
export PYTHONHASHSEED=0 PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONPYCACHEPREFIX=/dev/shm/orthohmm_native12_composed_driver_${SLURM_JOB_ID:?}
[[ ! -e "$PYTHONPYCACHEPREFIX" && ! -L "$PYTHONPYCACHEPREFIX" ]]
[[ "$#" == 4 && "$1" == --request && "$3" == --request-sha256 && "$4" == scheduler-comment ]]
REQUEST_SHA=$("$PYTHON" -B -c '
import re, subprocess, sys
from benchmark_tools.run_native12_composed_cost import SCRIPT, ROOT, JOB_NAME
from benchmark_tools.verify_threadripper_controller import validate
job = int(sys.argv[1])
result = subprocess.run(["scontrol", "show", "job", str(job), "--oneliner"],
    capture_output=True, text=True, check=True, timeout=5)
fields = validate(result.stdout, job, "running", command=str(SCRIPT), cwd=str(ROOT),
    time_limit="1-02:00:00", allocation_mode="shared")["fields"]
value = fields.get("Comment", "")
if fields.get("JobName") != JOB_NAME or not re.fullmatch("[0-9a-f]{64}", value):
    raise ValueError("Require exact native12 allocation and held request digest")
print(value)
' "$SLURM_JOB_ID")
exec "$PYTHON" -B benchmark_tools/run_native12_composed_cost.py --request "$2" --request-sha256 "$REQUEST_SHA"
