#!/usr/bin/env bash
set -euo pipefail

: "${TOOL:?TOOL must be orthofinder_v3_diamond or orthofinder_v3_mmseqs}"
: "${INPUT_DIR:?INPUT_DIR must contain the staged QfO FASTA files}"
: "${OUTDIR:?OUTDIR must name the result directory}"

case "${TOOL}" in
    orthofinder_v3_diamond) search=diamond ;;
    orthofinder_v3_mmseqs) search=mmseqs ;;
    *)
        echo "Unsupported tool: ${TOOL}" >&2
        exit 2
        ;;
esac

readonly SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
readonly BASE=$(dirname "${SCRIPT_DIR}")
readonly THREADS=${THREADS:-${SLURM_CPUS_PER_TASK:-32}}
readonly ANALYSIS_THREADS=${ANALYSIS_THREADS:-${THREADS}}
readonly PYTHON_BIN=${PYTHON_BIN:-python3}
readonly RUNDIR="${OUTDIR}/run"
readonly STAGE_DIR="${RUNDIR}/input"
readonly TIME_LOG="${OUTDIR}/time.log"
readonly QFO_PIPELINE=${QFO_PIPELINE:-${BASE}/qfo_benchmark/benchmark-webservice}
readonly QFO_MAPPING=${QFO_MAPPING:-${QFO_PIPELINE}/reference_data/2020/mapping.json.gz}

source "${SCRIPT_DIR}/orthofinder.sh"

if [[ $(find "${INPUT_DIR}" -maxdepth 1 -type f -name '*.fasta' | wc -l) -ne 78 ]]; then
    echo "Expected 78 QfO proteomes in ${INPUT_DIR}" >&2
    exit 1
fi

rm -rf "${RUNDIR}"
mkdir -p "${STAGE_DIR}"
cp -L "${INPUT_DIR}/"*.fasta "${STAGE_DIR}/"
orthofinder_verify_install > "${OUTDIR}/tool_version.txt"

start_epoch=$(date +%s)
if [[ "${search}" == "diamond" ]]; then
    /usr/bin/time -v -o "${TIME_LOG}" \
        "${ORTHOFINDER_BIN}" \
            -f "${STAGE_DIR}" -t "${THREADS}" -a "${ANALYSIS_THREADS}" -S diamond \
            > "${OUTDIR}/run.log" 2>&1
else
    readonly MMSEQS_BIN=${MMSEQS_BIN:-${SOFTWARE_ROOT}/mmseqs/bin/mmseqs}
    if [[ ! -x "${MMSEQS_BIN}" ]]; then
        echo "MMseqs executable not found: ${MMSEQS_BIN}" >&2
        exit 1
    fi
    /usr/bin/time -v -o "${TIME_LOG}" \
        env PATH="$(dirname "${MMSEQS_BIN}"):${PATH}" \
            OMP_NUM_THREADS="${THREADS}" MMSEQS_NUM_THREADS="${THREADS}" \
            "${ORTHOFINDER_BIN}" \
            -f "${STAGE_DIR}" -t "${THREADS}" -a "${ANALYSIS_THREADS}" -S mmseqs \
            > "${OUTDIR}/run.log" 2>&1
fi

results_dir=$(orthofinder_results_dir "${STAGE_DIR}")
orthofinder_collect_orthogroups "${results_dir}" "${OUTDIR}"
"${PYTHON_BIN}" "${SCRIPT_DIR}/orthofinder_to_pairwise.py" \
    "${results_dir}" > "${OUTDIR}/pairs.tsv" 2>> "${OUTDIR}/run.log"
"${PYTHON_BIN}" "${SCRIPT_DIR}/qfo_filter_pairs.py" \
    "${OUTDIR}/pairs.tsv" "${QFO_MAPPING}" "${OUTDIR}/pairs.qfo.tsv" \
    2>> "${OUTDIR}/run.log"

echo "$(( $(date +%s) - start_epoch ))" > "${OUTDIR}/wall_time_s.txt"
awk '/Maximum resident set size/ {print $NF}' "${TIME_LOG}" > "${OUTDIR}/peak_rss.txt"
