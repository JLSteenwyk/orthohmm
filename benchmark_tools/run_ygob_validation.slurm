#!/usr/bin/env bash
set -euo pipefail
[[ $# == 0 || ( $# == 1 && $1 == --check-only ) ]] || exit 2

: "${REPO_ROOT:?REPO_ROOT required}"
: "${FROZEN_ROOT:?FROZEN_ROOT required}"
: "${LAUNCHER_ROOT:?LAUNCHER_ROOT required}"
: "${LAUNCHER_COMMIT:?LAUNCHER_COMMIT required}"
: "${RUNTIME_MANIFEST:?RUNTIME_MANIFEST required}"
: "${RUNTIME_SHA256:?RUNTIME_SHA256 required}"
readonly SOURCE_COMMIT=7f3a9e40dd7e79f842cc2c11fb8b548f9a802806
readonly PREPARED="${REPO_ROOT}/benchmarks/work/ygob_validation_v1"
readonly OUT="${REPO_ROOT}/benchmarks/results/ygob_validation_v1"
readonly PYTHON=/home/bizon/anaconda3/bin/python
readonly OF=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/SOFTWARE/orthofinder_3.1.5/bin/orthofinder

[[ $(git -C "${FROZEN_ROOT}" rev-parse HEAD) == "${SOURCE_COMMIT}" ]]
git -C "${FROZEN_ROOT}" diff --quiet HEAD -- orthohmm benchmark_tools/benchmark_production.py
[[ $(git -C "${LAUNCHER_ROOT}" rev-parse HEAD) == "${LAUNCHER_COMMIT}" ]]
git -C "${LAUNCHER_ROOT}" diff --quiet HEAD -- benchmark_tools
cmp -s "${BASH_SOURCE[0]}" "${LAUNCHER_ROOT}/benchmark_tools/run_ygob_validation.slurm"
[[ ! -e "${OUT}" ]] || { echo "Output already exists: ${OUT}" >&2; exit 2; }
check_runtime() {
    [[ $(sha256sum "${RUNTIME_MANIFEST}" | cut -d ' ' -f 1) == "${RUNTIME_SHA256}" ]]
    "${PYTHON}" "${LAUNCHER_ROOT}/benchmark_tools/build_publication_runtime.py" \
        --root "${FROZEN_ROOT}" --output "${RUNTIME_MANIFEST}" --verify-only
}
check_runtime
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
if {p.resolve() for p in (root / "input").iterdir()} != {Path(r["path"]).resolve() for r in manifest["inputs"]}:
    raise SystemExit("Unexpected prepared input file set")
PY
if [[ ${1:-} == --check-only ]]; then
    "${PYTHON}" - "${LAUNCHER_ROOT}" "${FROZEN_ROOT}" <<'PY'
import sys
sys.path.insert(0, sys.argv[1])
from benchmark_tools.validate_profile_runtime import require_profile_runtime
require_profile_runtime(sys.argv[2])
print("YGOB input/source/native preflight passed; no inference")
PY
    exit 0
fi
mkdir -p "${OUT}"
trap 'printf "runner_exit_code\t%s\n" "$?" >> "${OUT}/run_metadata.tsv"' EXIT
printf 'job_id\t%s\nstarted\t%s\n' "${SLURM_JOB_ID:-local}" "$(date -u +%FT%TZ)" > "${OUT}/run_metadata.tsv"
git -C "${FROZEN_ROOT}" rev-parse HEAD > "${OUT}/inference_source_commit.txt"
cp -p "${BASH_SOURCE[0]}" "${OUT}/launcher.slurm"
cp -p "${RUNTIME_MANIFEST}" "${OUT}/native_runtime.json"
printf '%s\n' "${LAUNCHER_COMMIT}" > "${OUT}/launcher_commit.txt"
sha256sum "${PREPARED}/manifest.json" "${OUT}/launcher.slurm" "${OUT}/native_runtime.json" > "${OUT}/launcher_inputs.sha256"
"${PYTHON}" -c 'import importlib.metadata as m, json, sys; print(json.dumps({"python":sys.version,"packages":sorted((d.metadata["Name"],d.version) for d in m.distributions())},indent=2))' > "${OUT}/python_environment.json"
export PYTHONPATH="${FROZEN_ROOT}" PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PATH="/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/SOFTWARE/mafft-7.525-with-extensions/bin:/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/SOFTWARE/FastTree_v220:/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/SOFTWARE/diamond-linux64:${PATH}"
sha256sum "$(command -v mafft)" "$(command -v FastTree)" "$(command -v diamond)" "${OF}" > "${OUT}/tool_entrypoints.sha256"
"$(dirname "${OF}")/python" -c 'import importlib.metadata as m, json, sys; assert m.version("orthofinder") == "3.1.5"; print(json.dumps({"python":sys.version,"packages":sorted((d.metadata["Name"],d.version) for d in m.distributions())},indent=2))' > "${OUT}/orthofinder_environment.json"
runtime_probe() {
    check_runtime
    "${PYTHON}" "${LAUNCHER_ROOT}/benchmark_tools/validate_profile_runtime.py" \
        --root "${FROZEN_ROOT}" --output "${OUT}/${1}.runtime.json"
}

runtime_probe high_sensitivity_before
/usr/bin/time -v -o "${OUT}/high_sensitivity.time.log" \
    "${PYTHON}" "${FROZEN_ROOT}/benchmark_tools/benchmark_production.py" \
    "${PREPARED}/input" "${OUT}/high_sensitivity" "${OUT}/high_sensitivity.json" \
    --cpu 32 --threads-per-worker 8 --accuracy-profile high_sensitivity
runtime_probe high_sensitivity_after

runtime_probe satellite_v2_before
/usr/bin/time -v -o "${OUT}/satellite_v2.time.log" \
    "${PYTHON}" "${FROZEN_ROOT}/benchmark_tools/benchmark_production.py" \
    "${PREPARED}/input" "${OUT}/satellite_v2" "${OUT}/satellite_v2.json" \
    --cpu 32 --threads-per-worker 8 --accuracy-profile high_sensitivity \
    --phylogeny reconcile --species-tree-mode infer --phylogeny-candidates satellite_v2 \
    --phylogeny-root-rule species_overlap --phylogeny-pair-rule positive_paralogy \
    --species-tree-rooting min_variance
runtime_probe satellite_v2_after

mkdir -p "${OUT}/orthofinder/input"
cp -p "${PREPARED}/input/"*.fasta "${OUT}/orthofinder/input/"
sha256sum "${OF}" > "${OUT}/orthofinder_entrypoint.sha256"
/usr/bin/time -v -o "${OUT}/orthofinder.time.log" \
    "${OF}" -f "${OUT}/orthofinder/input" -t 32 -a 8 -S diamond \
    > "${OUT}/orthofinder.log" 2>&1
sha256sum --check --status "${OUT}/tool_entrypoints.sha256"
printf 'finished\t%s\n' "$(date -u +%FT%TZ)" >> "${OUT}/run_metadata.tsv"
# Accuracy scoring is deliberately a separate step gated on the overlap audit.
