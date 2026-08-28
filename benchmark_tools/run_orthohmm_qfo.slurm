#!/usr/bin/env bash
set -euo pipefail

: "${INPUT_DIR:?INPUT_DIR must contain the 78 QfO FASTA files}"
: "${OUTDIR:?OUTDIR must name a new result directory}"

: "${REPO_ROOT:?REPO_ROOT must point to the OrthoHMM repository}"

readonly BASE=${REPO_ROOT}
readonly SCRIPT_DIR="${BASE}/benchmark_tools"
readonly THREADS=${THREADS:-${SLURM_CPUS_PER_TASK:-32}}
readonly THREADS_PER_WORKER=${THREADS_PER_WORKER:-8}
readonly PYTHON_BIN=${PYTHON_BIN:-/home/bizon/anaconda3/bin/python}
readonly OUTPUT_DIR="${OUTDIR}/output"
readonly METRICS_JSON="${OUTDIR}/metrics.json"
readonly QFO_PIPELINE=${QFO_PIPELINE:-${BASE}/qfo_benchmark/benchmark-webservice}
readonly QFO_MAPPING=${QFO_MAPPING:-${QFO_PIPELINE}/reference_data/2020/mapping.json.gz}

if [[ $(find "${INPUT_DIR}" -maxdepth 1 -type f -name '*.fasta' | wc -l) -ne 78 ]]; then
    echo "Expected 78 QfO proteomes in ${INPUT_DIR}" >&2
    exit 1
fi
if [[ -e "${OUTPUT_DIR}" || -e "${METRICS_JSON}" ]]; then
    echo "Refusing to reuse an existing production result in ${OUTDIR}" >&2
    exit 1
fi

mkdir -p "${OUTDIR}"
git -C "${BASE}" rev-parse HEAD > "${OUTDIR}/source_commit.txt"
start_epoch=$(date +%s)

"${PYTHON_BIN}" "${SCRIPT_DIR}/benchmark_production.py" \
    "${INPUT_DIR}" \
    "${OUTPUT_DIR}" \
    "${METRICS_JSON}" \
    --cpu "${THREADS}" \
    --threads-per-worker "${THREADS_PER_WORKER}" \
    --accuracy-profile high_sensitivity \
    > "${OUTDIR}/run.log" 2>&1

cp "${OUTPUT_DIR}/orthohmm_orthogroups.txt" "${OUTDIR}/orthogroups.txt"
"${PYTHON_BIN}" "${BASE}/qfo_benchmark/og_to_pairwise.py" \
    "${OUTDIR}/orthogroups.txt" "${INPUT_DIR}" \
    > "${OUTDIR}/pairs.tsv" 2>> "${OUTDIR}/run.log"
"${PYTHON_BIN}" "${SCRIPT_DIR}/qfo_filter_pairs.py" \
    "${OUTDIR}/pairs.tsv" "${QFO_MAPPING}" "${OUTDIR}/pairs.qfo.tsv" \
    2>> "${OUTDIR}/run.log"

echo "$(( $(date +%s) - start_epoch ))" > "${OUTDIR}/wall_time_s.txt"
"${PYTHON_BIN}" - "${METRICS_JSON}" "${OUTDIR}/peak_rss.txt" <<'PY'
import json
import sys

metrics = json.load(open(sys.argv[1]))
peak_kib = int(metrics["peak_process_tree_rss_bytes"]) // 1024
with open(sys.argv[2], "w") as handle:
    handle.write(f"{peak_kib}\n")
PY
