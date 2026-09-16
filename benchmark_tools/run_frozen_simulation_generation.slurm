#!/usr/bin/env bash
#SBATCH --job-name=variable_sim_generation
#SBATCH --partition=gpu
#SBATCH --nodelist=bizon
#SBATCH --array=0-39%2
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=12:00:00
set -euo pipefail

: "${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID required}"
readonly EXECUTOR="${1:?Pinned worktree required}"
readonly COMMIT="${2:?Pinned commit required}"
readonly MANIFEST="${3:?Frozen generation manifest required}"
readonly HASH="${4:?Manifest SHA256 required}"
readonly PANEL="${5:?Panel root required}"
readonly PYTHON=/home/bizon/anaconda3/bin/python
shift 5
[[ $# == 0 || ( $# == 1 && $1 == --check-only ) ]] || exit 2
[[ $(git -C "${EXECUTOR}" rev-parse HEAD) == "${COMMIT}" ]]
git -C "${EXECUTOR}" diff --quiet HEAD -- benchmark_tools
label=$("${PYTHON}" - "${MANIFEST}" "${HASH}" "${SLURM_ARRAY_TASK_ID}" <<'PY'
import hashlib
import json
from pathlib import Path
import sys
raw = Path(sys.argv[1]).read_bytes()
if hashlib.sha256(raw).hexdigest() != sys.argv[2]:
    raise SystemExit("Generation manifest changed")
runs = json.loads(raw)["simulation_runs"]
index = int(sys.argv[3])
if len(runs) != 40 or not 0 <= index < 40:
    raise SystemExit("Invalid generation array index")
print(runs[index]["label"])
PY
)
printf 'array_job=%s task=%s run=%s executor_commit=%s\n' \
    "${SLURM_ARRAY_JOB_ID:-check}" "${SLURM_ARRAY_TASK_ID}" "${label}" "${COMMIT}"
exec "${PYTHON}" "${EXECUTOR}/benchmark_tools/run_simulation_generation.py" \
    --manifest "${MANIFEST}" --manifest-sha256 "${HASH}" --panel "${PANEL}" --run "${label}" \
    --evidence "${PANEL}/execution/${label}" "$@"
