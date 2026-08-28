#!/usr/bin/env bash

# Shared OrthoFinder configuration for every benchmark runner.
readonly ORTHOFINDER_VERSION="3.1.5"
readonly ORTHOFINDER_SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
readonly ORTHOFINDER_REPO_ROOT=$(dirname "${ORTHOFINDER_SCRIPT_DIR}")
readonly ORTHOFINDER_WORKSPACE_ROOT=$(dirname "$(dirname "${ORTHOFINDER_REPO_ROOT}")")
SOFTWARE_ROOT=${SOFTWARE_ROOT:-${ORTHOFINDER_WORKSPACE_ROOT}/SOFTWARE}
readonly ORTHOFINDER_DEFAULT_BIN="${SOFTWARE_ROOT}/orthofinder_3.1.5/bin/orthofinder"
ORTHOFINDER_BIN="${ORTHOFINDER_BIN:-${ORTHOFINDER_DEFAULT_BIN}}"

orthofinder_verify_install() {
    if [[ ! -x "${ORTHOFINDER_BIN}" ]]; then
        echo "OrthoFinder ${ORTHOFINDER_VERSION} executable not found: ${ORTHOFINDER_BIN}" >&2
        return 1
    fi

    local detected
    detected=$("${ORTHOFINDER_BIN}" --version 2>&1)
    if [[ "${detected}" != "OrthoFinder:v${ORTHOFINDER_VERSION}" ]]; then
        echo "Expected OrthoFinder:v${ORTHOFINDER_VERSION}, got: ${detected}" >&2
        return 1
    fi
    printf '%s\n' "${detected}"
}

orthofinder_results_dir() {
    local input_dir=$1
    local candidates=()
    mapfile -t candidates < <(
        find "${input_dir}/OrthoFinder" -type f \
            -path '*/Orthogroups/Orthogroups.txt' -print 2>/dev/null
    )
    if [[ ${#candidates[@]} -ne 1 ]]; then
        echo "Expected one OrthoFinder result under ${input_dir}; found ${#candidates[@]}" >&2
        return 1
    fi
    dirname "$(dirname "${candidates[0]}")"
}

orthofinder_collect_orthogroups() {
    local results_dir=$1
    local output_dir=$2
    local log_file="${results_dir}/Log.txt"
    local source_file="${results_dir}/Orthogroups/Orthogroups.txt"

    if ! grep -Fq "Started OrthoFinder version ${ORTHOFINDER_VERSION}" "${log_file}"; then
        echo "Result provenance does not report OrthoFinder ${ORTHOFINDER_VERSION}: ${log_file}" >&2
        return 1
    fi
    if [[ ! -s "${source_file}" ]]; then
        echo "Missing OrthoFinder root orthogroups: ${source_file}" >&2
        return 1
    fi

    cp "${source_file}" "${output_dir}/orthogroups.txt"
    printf 'OrthoFinder:v%s\n' "${ORTHOFINDER_VERSION}" > "${output_dir}/tool_version.txt"
    printf '%s\n' \
        'Orthogroups/Orthogroups.txt (v3.1.5 root HOGs; equivalent to pre-3.1.5 N0.tsv)' \
        > "${output_dir}/orthogroups_source.txt"
}
