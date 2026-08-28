#!/usr/bin/env bash
set -euo pipefail

readonly VERSION=3.1.5
readonly ARCHIVE_NAME="orthofinder-linux-intel-${VERSION}.tar.gz"
readonly ARCHIVE_SHA256=d74b5dbf9348e7ffdb735f1fda4ff9c49f504261811355519cf9dafcf80e7739
readonly DOWNLOAD_URL="https://github.com/OrthoFinder/OrthoFinder/releases/download/v${VERSION}/${ARCHIVE_NAME}"
readonly SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
readonly REPO_ROOT=$(dirname "${SCRIPT_DIR}")
readonly WORKSPACE_ROOT=$(dirname "$(dirname "${REPO_ROOT}")")
readonly SOFTWARE_DIR=${SOFTWARE_ROOT:-${WORKSPACE_ROOT}/SOFTWARE}
readonly INSTALL_DIR=${1:-${SOFTWARE_DIR}/orthofinder_${VERSION}}
readonly ARCHIVE_DIR="${SOFTWARE_DIR}/downloads"
readonly ARCHIVE="${ARCHIVE_DIR}/${ARCHIVE_NAME}"

mkdir -p "${ARCHIVE_DIR}"
if [[ ! -f "${ARCHIVE}" ]] || ! echo "${ARCHIVE_SHA256}  ${ARCHIVE}" | sha256sum -c -; then
    curl -fL --retry 3 -o "${ARCHIVE}" "${DOWNLOAD_URL}"
fi
echo "${ARCHIVE_SHA256}  ${ARCHIVE}" | sha256sum -c -

rm -rf "${INSTALL_DIR}"
python3.12 -m venv "${INSTALL_DIR}"
"${INSTALL_DIR}/bin/python" -m pip install --disable-pip-version-check "${ARCHIVE}"

detected=$("${INSTALL_DIR}/bin/orthofinder" --version)
if [[ "${detected}" != "OrthoFinder:v${VERSION}" ]]; then
    echo "Installed executable reports an unexpected version: ${detected}" >&2
    exit 1
fi

"${INSTALL_DIR}/bin/python" -m pip freeze | sort > "${INSTALL_DIR}/requirements.lock"
printf '%s  %s\n' "${ARCHIVE_SHA256}" "${ARCHIVE_NAME}" > "${INSTALL_DIR}/source.sha256"
printf 'Installed %s at %s\n' "${detected}" "${INSTALL_DIR}"
