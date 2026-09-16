#!/usr/bin/env bash
#SBATCH --job-name=simulation_generation
#SBATCH --partition=gpu
#SBATCH --nodelist=bizon
#SBATCH --array=0-39%2
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=12:00:00
set -euo pipefail

: "${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID required}"
readonly ROOT="${1:?Repository root required}"
readonly PYTHON=/home/bizon/anaconda3/bin/python
readonly MANIFEST="${ROOT}/benchmark_tools/results/publication_simulation_manifest_20260916.json"
readonly MANIFEST_SHA=ee31ea38d3b5c047abf80f04636f06838959f941d6704649a216bb93958343b2
readonly RUNNER="${ROOT}/benchmarks/work/publication_methods_frozen/benchmark_tools/run_simulation_generation.py"
readonly PANEL="${ROOT}/benchmarks/work/publication_simulation_panel_v1"
shift
[[ $# == 0 || ( $# == 1 && $1 == --check-only ) ]] || exit 2

label=$("${PYTHON}" - "${MANIFEST}" "${MANIFEST_SHA}" "${RUNNER}" "${SLURM_ARRAY_TASK_ID}" <<'PY'
import hashlib
import json
from pathlib import Path
import sys
manifest, expected, runner, index = sys.argv[1:]
raw = Path(manifest).read_bytes()
if hashlib.sha256(raw).hexdigest() != expected:
    raise SystemExit("Generation manifest changed")
if hashlib.sha256(Path(runner).read_bytes()).hexdigest() != "00eb976fb5a3de1ec9badc196dcfc48c8bcf110beb447eb6f878f7fb6ae90f91":
    raise SystemExit("Pinned generation runner changed")
runs = json.loads(raw)["simulation_runs"]
if len(runs) != 40 or not 0 <= int(index) < len(runs):
    raise SystemExit("Invalid generation array index")
print(runs[int(index)]["label"])
PY
)
printf 'array_job=%s task=%s run=%s\n' "${SLURM_ARRAY_JOB_ID:-check}" "${SLURM_ARRAY_TASK_ID}" "${label}"
exec "${PYTHON}" "${RUNNER}" --manifest "${MANIFEST}" \
    --manifest-sha256 "${MANIFEST_SHA}" --panel "${PANEL}" --run "${label}" \
    --evidence "${PANEL}/execution/${label}" "$@"
