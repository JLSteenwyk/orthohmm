#!/usr/bin/env bash
#SBATCH --job-name=simulation_methods
#SBATCH --partition=gpu
#SBATCH --nodelist=bizon
#SBATCH --array=0-69%2
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
set -euo pipefail

: "${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID required}"
readonly ROOT="${1:?Repository root required}"
readonly EXECUTOR="${2:?Pinned executor worktree required}"
readonly COMMIT="${3:?Pinned executor commit required}"
readonly PYTHON=/home/bizon/anaconda3/bin/python
readonly MANIFEST="${ROOT}/benchmark_tools/results/publication_simulation_methods_20260916.json"
shift 3
[[ $# == 0 || ( $# == 1 && $1 == --check-only ) ]] || exit 2
[[ $(git -C "${EXECUTOR}" rev-parse HEAD) == "${COMMIT}" ]]
git -C "${EXECUTOR}" diff --quiet HEAD -- benchmark_tools
label=$("${PYTHON}" - "${MANIFEST}" "${SLURM_ARRAY_TASK_ID}" <<'PY'
import hashlib
import json
from pathlib import Path
import sys
raw = Path(sys.argv[1]).read_bytes()
if hashlib.sha256(raw).hexdigest() != "4f0717fa8b0c34d5a6982a3193772d64bc39a3d75b8fe23195d9ce1b21eee186":
    raise SystemExit("Method manifest changed")
datasets = json.loads(raw)["datasets"]
index = int(sys.argv[2])
if len(datasets) != 70 or not 0 <= index < 70:
    raise SystemExit("Invalid method array index")
print(datasets[index]["label"])
PY
)
printf 'array_job=%s task=%s dataset=%s executor_commit=%s\n' \
    "${SLURM_ARRAY_JOB_ID:-check}" "${SLURM_ARRAY_TASK_ID}" "${label}" "${COMMIT}"
exec "${PYTHON}" "${EXECUTOR}/benchmark_tools/run_simulation_methods.py" \
    --manifest "${MANIFEST}" --panel "${ROOT}/benchmarks/work/publication_simulation_panel_v1" \
    --dataset "${label}" --evidence "${ROOT}/benchmarks/results/publication_simulation_methods_v1/${label}/execution" "$@"
