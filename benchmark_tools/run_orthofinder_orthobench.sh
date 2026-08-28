#!/usr/bin/env bash
set -euo pipefail

readonly SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
readonly BASE=$(dirname "${SCRIPT_DIR}")
readonly WORKSPACE_ROOT=$(dirname "$(dirname "${BASE}")")
readonly ORTHOBENCH=${ORTHOBENCH_ROOT:-${WORKSPACE_ROOT}/Open_Orthobench/BENCHMARKS}
readonly SEARCH=${1:-diamond}
readonly THREADS=${THREADS:-32}
readonly ANALYSIS_THREADS=${ANALYSIS_THREADS:-8}
readonly OUTDIR="${BASE}/benchmarks/results/orthofinder_v3_${SEARCH}"
readonly STAGE_DIR="${OUTDIR}/input"

source "${BASE}/benchmark_tools/orthofinder.sh"

if [[ "${SEARCH}" != "diamond" && "${SEARCH}" != "mmseqs" ]]; then
    echo "Usage: $0 [diamond|mmseqs]" >&2
    exit 2
fi

rm -rf "${OUTDIR}"
mkdir -p "${STAGE_DIR}"
cp "${ORTHOBENCH}/Input/"*.fa "${STAGE_DIR}/"
orthofinder_verify_install > "${OUTDIR}/tool_version.txt"

if [[ "${SEARCH}" == "diamond" ]]; then
    /usr/bin/time -v -o "${OUTDIR}/time.log" \
        "${ORTHOFINDER_BIN}" \
            -f "${STAGE_DIR}" -t "${THREADS}" -a "${ANALYSIS_THREADS}" -S diamond \
            > "${OUTDIR}/run.log" 2>&1
else
    readonly MMSEQS_BIN=${MMSEQS_BIN:-${SOFTWARE_ROOT}/mmseqs/bin/mmseqs}
    if [[ ! -x "${MMSEQS_BIN}" ]]; then
        echo "MMseqs executable not found: ${MMSEQS_BIN}" >&2
        exit 1
    fi
    /usr/bin/time -v -o "${OUTDIR}/time.log" \
        env \
            PATH="$(dirname "${MMSEQS_BIN}"):${PATH}" \
            OMP_NUM_THREADS="${THREADS}" \
            MMSEQS_NUM_THREADS="${THREADS}" \
            "${ORTHOFINDER_BIN}" \
            -f "${STAGE_DIR}" -t "${THREADS}" -a "${ANALYSIS_THREADS}" -S mmseqs \
            > "${OUTDIR}/run.log" 2>&1
fi

OF_RESULTS=$(orthofinder_results_dir "${STAGE_DIR}")
orthofinder_collect_orthogroups "${OF_RESULTS}" "${OUTDIR}"
cp "${OUTDIR}/orthogroups.txt" \
    "${BASE}/benchmarks/results/orthogroups_orthofinder_v3_${SEARCH}.txt"

python "${ORTHOBENCH}/benchmark.py" "${OUTDIR}/orthogroups.txt" \
    | tee "${OUTDIR}/orthobench_score.txt"
