#!/usr/bin/env bash
set -euo pipefail

: "${REPO_ROOT:?REPO_ROOT required}"
: "${FROZEN_ROOT:?FROZEN_ROOT required}"
readonly SOURCE_COMMIT=7f3a9e4
readonly PREPARED="${REPO_ROOT}/benchmarks/work/ygob_validation_v1"
readonly OUT="${REPO_ROOT}/benchmarks/results/ygob_validation_v1"
readonly PYTHON=/home/bizon/anaconda3/bin/python
readonly OF=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/SOFTWARE/orthofinder_3.1.5/bin/orthofinder

[[ $(git -C "${FROZEN_ROOT}" rev-parse --short=7 HEAD) == "${SOURCE_COMMIT}" ]]
git -C "${FROZEN_ROOT}" diff --quiet HEAD -- orthohmm benchmark_tools/benchmark_production.py
[[ ! -e "${OUT}" ]] || { echo "Output already exists: ${OUT}" >&2; exit 2; }
"${PYTHON}" - "${PREPARED}" <<'PY'
import hashlib
import json
from pathlib import Path
import sys
root = Path(sys.argv[1])
manifest = json.loads((root / "manifest.json").read_text())
if manifest["proteins"] != 83404 or len(manifest["species"]) != 16:
    raise SystemExit("Unexpected validation dimensions")
for item in manifest["inputs"] + [manifest["reference"]]:
    path = Path(item["path"])
    if hashlib.sha256(path.read_bytes()).hexdigest() != item["sha256"]:
        raise SystemExit(f"Changed prepared input: {path}")
PY
mkdir -p "${OUT}"
trap 'printf "runner_exit_code\t%s\n" "$?" >> "${OUT}/run_metadata.tsv"' EXIT
printf 'job_id\t%s\nstarted\t%s\n' "${SLURM_JOB_ID:-local}" "$(date -u +%FT%TZ)" > "${OUT}/run_metadata.tsv"
git -C "${FROZEN_ROOT}" rev-parse HEAD > "${OUT}/inference_source_commit.txt"
sha256sum "${PREPARED}/manifest.json" "${BASH_SOURCE[0]}" > "${OUT}/launcher_inputs.sha256"
"${PYTHON}" -c 'import importlib.metadata as m, json, sys; print(json.dumps({"python":sys.version,"packages":sorted((d.metadata["Name"],d.version) for d in m.distributions())},indent=2))' > "${OUT}/python_environment.json"
export PYTHONPATH="${FROZEN_ROOT}" PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

/usr/bin/time -v -o "${OUT}/high_sensitivity.time.log" \
    "${PYTHON}" "${FROZEN_ROOT}/benchmark_tools/benchmark_production.py" \
    "${PREPARED}/input" "${OUT}/high_sensitivity" "${OUT}/high_sensitivity.json" \
    --cpu 32 --threads-per-worker 8 --accuracy-profile high_sensitivity

/usr/bin/time -v -o "${OUT}/satellite_v2.time.log" \
    "${PYTHON}" "${FROZEN_ROOT}/benchmark_tools/benchmark_production.py" \
    "${PREPARED}/input" "${OUT}/satellite_v2" "${OUT}/satellite_v2.json" \
    --cpu 32 --threads-per-worker 8 --accuracy-profile high_sensitivity \
    --phylogeny reconcile --species-tree-mode infer --phylogeny-candidates satellite_v2 \
    --phylogeny-root-rule species_overlap --phylogeny-pair-rule positive_paralogy \
    --species-tree-rooting min_variance

mkdir -p "${OUT}/orthofinder/input"
cp -p "${PREPARED}/input/"*.fasta "${OUT}/orthofinder/input/"
sha256sum "${OF}" > "${OUT}/orthofinder_entrypoint.sha256"
/usr/bin/time -v -o "${OUT}/orthofinder.time.log" \
    "${OF}" -f "${OUT}/orthofinder/input" -t 32 -a 8 -S diamond \
    > "${OUT}/orthofinder.log" 2>&1
printf 'finished\t%s\n' "$(date -u +%FT%TZ)" >> "${OUT}/run_metadata.tsv"
# Accuracy scoring is deliberately a separate step gated on the overlap audit.
