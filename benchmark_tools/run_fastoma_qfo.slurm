#!/usr/bin/env bash
set -euo pipefail

: "${REPO_ROOT:?REPO_ROOT must point to the OrthoHMM repository}"
: "${INPUT_DIR:?INPUT_DIR must contain the 78 QfO FASTA files}"
: "${OUTDIR:?OUTDIR must name a new result directory}"

readonly BASE=${REPO_ROOT}
readonly SCRIPT_DIR="${BASE}/benchmark_tools"
readonly FASTOMA_ROOT=${FASTOMA_ROOT:-/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/SOFTWARE/FastOMA}
readonly FASTOMA_VERSION=${FASTOMA_VERSION:-0.3.5}
readonly NEXTFLOW_BIN=${NEXTFLOW_BIN:-/home/bizon/bin/nextflow}
readonly JAVA_HOME_17=${JAVA_HOME_17:-/home/bizon/anaconda3/pkgs/openjdk-17.0.18-ha668962_0/lib/jvm}
readonly THREADS=${THREADS:-${SLURM_CPUS_PER_TASK:-32}}
readonly MAX_MEMORY=${MAX_MEMORY:-700.GB}
readonly RUN_DIR="${OUTDIR}/run"
readonly STAGE_DIR="${RUN_DIR}/input"
readonly PROTEOME_DIR="${STAGE_DIR}/proteome"
readonly OUTPUT_DIR="${OUTDIR}/output"
readonly WORK_DIR="${OUTDIR}/work"
readonly TIME_LOG="${OUTDIR}/time.log"
readonly QFO_PIPELINE=${QFO_PIPELINE:-${BASE}/qfo_benchmark/benchmark-webservice}
readonly QFO_MAPPING=${QFO_MAPPING:-${QFO_PIPELINE}/reference_data/2020/mapping.json.gz}
readonly SPECIES_TREE=${SPECIES_TREE:-${BASE}/qfo_benchmark/results/orthofinder_v3_diamond/run/input/OrthoFinder/Results_Jul15/Species_Tree/SpeciesTree_rooted.txt}
readonly OMAMER_DB=${OMAMER_DB:-${BASE}/qfo_benchmark/data/LUCA.h5}

if [[ $(find "${INPUT_DIR}" -maxdepth 1 -type f -name '*.fasta' | wc -l) -ne 78 ]]; then
    echo "Expected 78 QfO proteomes in ${INPUT_DIR}" >&2
    exit 1
fi
for required in "${NEXTFLOW_BIN}" "${FASTOMA_ROOT}/FastOMA.nf" "${SPECIES_TREE}" "${OMAMER_DB}" "${QFO_MAPPING}"; do
    if [[ ! -s "${required}" ]]; then
        echo "Required FastOMA input is missing or empty: ${required}" >&2
        exit 1
    fi
done
if [[ -e "${OUTPUT_DIR}" || -e "${OUTDIR}/pairs.qfo.tsv" ]]; then
    echo "Refusing to reuse an existing production result in ${OUTDIR}" >&2
    exit 1
fi

mkdir -p "${PROTEOME_DIR}" "${OUTDIR}"
for source in "${INPUT_DIR}/"*.fasta; do
    filename=$(basename "${source}" .fasta)
    cp -L --reflink=auto "${source}" "${PROTEOME_DIR}/${filename}.fa"
done
cp -L "${SPECIES_TREE}" "${STAGE_DIR}/species_tree.nwk"

printf 'FastOMA version: %s\n' "${FASTOMA_VERSION}" > "${OUTDIR}/tool_version.txt"
printf 'FastOMA source: %s\n' "$(git -C "${FASTOMA_ROOT}" rev-parse HEAD)" >> "${OUTDIR}/tool_version.txt"
sha256sum "${STAGE_DIR}/species_tree.nwk" > "${OUTDIR}/species_tree.sha256"

start_epoch=$(date +%s)
cd "${FASTOMA_ROOT}"
export JAVA_HOME="${JAVA_HOME_17}"
export JAVA_CMD="${JAVA_HOME_17}/bin/java"
/usr/bin/time -v -o "${TIME_LOG}" \
    "${NEXTFLOW_BIN}" run FastOMA.nf \
        -profile docker \
        -ansi-log false \
        --container_version "${FASTOMA_VERSION}" \
        --input_folder "${STAGE_DIR}" \
        --output_folder "${OUTPUT_DIR}" \
        --omamer_db "${OMAMER_DB}" \
        --fasta_header_id_transformer uniprot \
        --force_pairwise_ortholog_generation true \
        --max_cpus "${THREADS}" \
        --max_memory "${MAX_MEMORY}" \
        -work-dir "${WORK_DIR}" \
        > "${OUTDIR}/run.log" 2>&1

if [[ ! -s "${OUTPUT_DIR}/orthologs.tsv.gz" ]]; then
    echo "FastOMA did not produce pairwise ortholog predictions" >&2
    exit 1
fi

gzip -cd "${OUTPUT_DIR}/orthologs.tsv.gz" \
    | awk -F '\t' 'BEGIN {OFS="\t"} NF == 2 {if ($1 < $2) print $1, $2; else print $2, $1}' \
    > "${OUTDIR}/pairs.tsv"
python3 "${SCRIPT_DIR}/qfo_filter_pairs.py" \
    "${OUTDIR}/pairs.tsv" "${QFO_MAPPING}" "${OUTDIR}/pairs.qfo.tsv" \
    2>> "${OUTDIR}/run.log"

echo "$(( $(date +%s) - start_epoch ))" > "${OUTDIR}/wall_time_s.txt"
awk '/Maximum resident set size/ {print $NF}' "${TIME_LOG}" > "${OUTDIR}/peak_rss.txt"
wc -l "${OUTDIR}/pairs.tsv" "${OUTDIR}/pairs.qfo.tsv" > "${OUTDIR}/pair_counts.txt"
