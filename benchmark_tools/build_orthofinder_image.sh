#!/usr/bin/env bash
set -euo pipefail

readonly SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
readonly REPO_ROOT=$(dirname "${SCRIPT_DIR}")
readonly WORKSPACE_ROOT=$(dirname "$(dirname "${REPO_ROOT}")")
readonly SOFTWARE_DIR=${SOFTWARE_ROOT:-${WORKSPACE_ROOT}/SOFTWARE}
readonly OUTPUT=${1:-${SOFTWARE_DIR}/orthofinder_3.1.5.sif}

mkdir -p "$(dirname "${OUTPUT}")"
if [[ ${EUID} -ne 0 ]]; then
    echo "This Singularity installation requires a privileged local build." >&2
    echo "Run this script as root, or use benchmark_tools/install_orthofinder.sh instead." >&2
    exit 1
fi
singularity build "${OUTPUT}" "${SCRIPT_DIR}/orthofinder-3.1.5.def"

detected=$(singularity exec "${OUTPUT}" orthofinder --version)
if [[ "${detected}" != "OrthoFinder:v3.1.5" ]]; then
    echo "Built image reports an unexpected version: ${detected}" >&2
    exit 1
fi

sha256sum "${OUTPUT}" > "${OUTPUT}.sha256"
printf 'Built %s (%s)\n' "${OUTPUT}" "${detected}"
