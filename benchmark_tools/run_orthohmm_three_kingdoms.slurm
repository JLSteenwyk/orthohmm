#!/usr/bin/env bash
set -euo pipefail

: "${REPO_ROOT:?REPO_ROOT must point to the OrthoHMM repository}"
: "${INPUT_DIR:?INPUT_DIR must contain the 12 Three Kingdoms FASTA files}"
: "${REFERENCE:?REFERENCE must point to the BUSCO reference orthogroups}"
: "${OUTDIR:?OUTDIR must name a new result directory}"

readonly BASE=${REPO_ROOT}
readonly THREADS=${THREADS:-${SLURM_CPUS_PER_TASK:-32}}
readonly THREADS_PER_WORKER=${THREADS_PER_WORKER:-8}
readonly PYTHON_BIN=${PYTHON_BIN:-/home/bizon/anaconda3/bin/python}
readonly OUTPUT_DIR="${OUTDIR}/output"
readonly METRICS_JSON="${OUTDIR}/metrics.json"
readonly SCORER=${SCORER:-$(dirname "${INPUT_DIR}")/score_against_busco.py}

if [[ $(find "${INPUT_DIR}" -maxdepth 1 -type f -name '*.fasta' | wc -l) -ne 12 ]]; then
    echo "Expected 12 Three Kingdoms proteomes in ${INPUT_DIR}" >&2
    exit 1
fi
if [[ ! -f "${REFERENCE}" ]]; then
    echo "BUSCO reference not found: ${REFERENCE}" >&2
    exit 1
fi
if [[ ! -f "${SCORER}" ]]; then
    echo "Three Kingdoms scorer not found: ${SCORER}" >&2
    exit 1
fi
if [[ -e "${OUTPUT_DIR}" || -e "${METRICS_JSON}" ]]; then
    echo "Refusing to reuse an existing production result in ${OUTDIR}" >&2
    exit 1
fi

mkdir -p "${OUTDIR}"
git -C "${BASE}" rev-parse HEAD > "${OUTDIR}/source_commit.txt"

"${PYTHON_BIN}" "${BASE}/benchmark_tools/benchmark_production.py" \
    "${INPUT_DIR}" \
    "${OUTPUT_DIR}" \
    "${METRICS_JSON}" \
    --cpu "${THREADS}" \
    --threads-per-worker "${THREADS_PER_WORKER}" \
    --accuracy-profile high_sensitivity \
    > "${OUTDIR}/run.log" 2>&1

cp "${OUTPUT_DIR}/orthohmm_orthogroups.txt" "${OUTDIR}/orthogroups.txt"
"${PYTHON_BIN}" "${SCORER}" \
    --predictions "${OUTDIR}/orthogroups.txt" \
    --reference "${REFERENCE}" \
    --label orthohmm_high_sensitivity \
    > "${OUTDIR}/score.txt"
