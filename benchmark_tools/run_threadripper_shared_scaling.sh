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
if [[ "$#" == 4 && "$1" == --request && "$3" == --request-sha256 && "$4" == scheduler-comment ]]; then
    REQUEST_SHA=$("$PYTHON" -B -c '
import re
import subprocess
import sys
from benchmark_tools.verify_threadripper_controller import validate
job = int(sys.argv[1])
result = subprocess.run(["scontrol", "show", "job", str(job), "--oneliner"],
    capture_output=True, text=True, check=True, timeout=5)
fields = validate(result.stdout, job, "running", command=sys.argv[2], cwd=sys.argv[3],
    time_limit="1-02:00:00", allocation_mode="shared")["fields"]
value = fields.get("Comment", "")
if not re.fullmatch("[0-9a-f]{64}", value):
    raise ValueError("Require the actual held-job request hash in its scheduler comment")
print(value)
' "$SLURM_JOB_ID" "$ROOT/benchmark_tools/run_threadripper_shared_scaling.sh" "$ROOT")
    set -- "$1" "$2" "$3" "$REQUEST_SHA"
fi
exec "$PYTHON" -B benchmark_tools/run_threadripper_scaling.py "$@"
