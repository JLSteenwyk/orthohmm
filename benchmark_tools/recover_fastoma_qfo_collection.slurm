#!/usr/bin/env bash
set -euo pipefail

: "${REPO_ROOT:?REPO_ROOT must point to the OrthoHMM repository}"
: "${OUTDIR:?OUTDIR must point to the existing FastOMA result directory}"
: "${COLLECT_WORK_DIR:?COLLECT_WORK_DIR must point to a staged collect_subhogs work directory}"

readonly SCRIPT_DIR="${REPO_ROOT}/benchmark_tools"
readonly WORK_DIR="${OUTDIR}/work"
readonly OUTPUT_DIR="${OUTDIR}/output"
readonly QFO_PIPELINE=${QFO_PIPELINE:-${REPO_ROOT}/qfo_benchmark/benchmark-webservice}
readonly QFO_MAPPING=${QFO_MAPPING:-${QFO_PIPELINE}/reference_data/2020/mapping.json.gz}
readonly FASTOMA_VERSION=${FASTOMA_VERSION:-0.3.5}
readonly FASTOMA_IMAGE=${FASTOMA_IMAGE:-dessimozlab/fastoma:${FASTOMA_VERSION}}
readonly THREADS=${THREADS:-${SLURM_CPUS_PER_TASK:-16}}
readonly DOCKER_MEMORY=${DOCKER_MEMORY:-280g}

for required in \
    "${COLLECT_WORK_DIR}/pickle_folders" \
    "${COLLECT_WORK_DIR}/omamer_rhogs" \
    "${COLLECT_WORK_DIR}/gene_id_dic_xml.pickle" \
    "${COLLECT_WORK_DIR}/species_tree_checked.nwk" \
    "${QFO_MAPPING}"; do
    if [[ ! -e "${required}" ]]; then
        echo "Required recovery input is missing: ${required}" >&2
        exit 1
    fi
done
if [[ $(find "${COLLECT_WORK_DIR}/pickle_folders" -mindepth 1 -maxdepth 1 | wc -l) -ne 271 ]]; then
    echo "Expected all 271 FastOMA sub-HOG pickle batches" >&2
    exit 1
fi
if find "${COLLECT_WORK_DIR}" -xtype l -print -quit | grep -q .; then
    echo "The staged FastOMA collection inputs contain broken links" >&2
    exit 1
fi
if [[ -e "${OUTPUT_DIR}/orthologs.tsv.gz" || -e "${OUTDIR}/pairs.qfo.tsv" ]]; then
    echo "Refusing to overwrite an existing FastOMA pair result" >&2
    exit 1
fi

for output in FastOMA_HOGs.orthoxml OrthologousGroupsFasta OrthologousGroups.tsv \
    RootHOGs.tsv RootHOGsFasta orthologs.tsv.gz; do
    if [[ -e "${COLLECT_WORK_DIR}/${output}" ]]; then
        echo "Refusing to overwrite partial recovery output: ${COLLECT_WORK_DIR}/${output}" >&2
        exit 1
    fi
done

docker run --rm \
    --cpus "${THREADS}" \
    --memory "${DOCKER_MEMORY}" \
    -v "${WORK_DIR}:${WORK_DIR}" \
    -w "${COLLECT_WORK_DIR}" \
    "${FASTOMA_IMAGE}" \
    fastoma-collect-subhogs \
        --pickle-folder pickle_folders/ \
        --roothogs-folder omamer_rhogs/ \
        --gene-id-pickle-file gene_id_dic_xml.pickle \
        --out FastOMA_HOGs.orthoxml \
        --marker-groups-fasta OrthologousGroups.tsv \
        --roothog-tsv RootHOGs.tsv \
        --species-tree species_tree_checked.nwk \
        --id-transform UniProt \
        -vv \
    > "${OUTDIR}/collection_recovery.log" 2>&1

docker run --rm \
    --cpus "${THREADS}" \
    --memory "${DOCKER_MEMORY}" \
    -v "${WORK_DIR}:${WORK_DIR}" \
    -w "${COLLECT_WORK_DIR}" \
    "${FASTOMA_IMAGE}" \
    fastoma-helper -vv pw-rel \
        --orthoxml FastOMA_HOGs.orthoxml \
        --out orthologs.tsv.gz \
        --type ortholog \
    >> "${OUTDIR}/collection_recovery.log" 2>&1

mkdir -p "${OUTPUT_DIR}"
for output in FastOMA_HOGs.orthoxml OrthologousGroupsFasta OrthologousGroups.tsv \
    RootHOGs.tsv RootHOGsFasta orthologs.tsv.gz; do
    cp -a "${COLLECT_WORK_DIR}/${output}" "${OUTPUT_DIR}/"
done
gzip -t "${OUTPUT_DIR}/orthologs.tsv.gz"

gzip -cd "${OUTPUT_DIR}/orthologs.tsv.gz" \
    | awk -F '\t' 'BEGIN {OFS="\t"} NF == 2 {if ($1 < $2) print $1, $2; else print $2, $1}' \
    > "${OUTDIR}/pairs.tsv"
python3 "${SCRIPT_DIR}/qfo_filter_pairs.py" \
    "${OUTDIR}/pairs.tsv" "${QFO_MAPPING}" "${OUTDIR}/pairs.qfo.tsv" \
    2>> "${OUTDIR}/collection_recovery.log"

if [[ -s "${OUTDIR}/run_start_epoch.txt" ]]; then
    start_epoch=$(<"${OUTDIR}/run_start_epoch.txt")
    echo "$(( $(date +%s) - start_epoch ))" > "${OUTDIR}/wall_time_s.txt"
fi
wc -l "${OUTDIR}/pairs.tsv" "${OUTDIR}/pairs.qfo.tsv" > "${OUTDIR}/pair_counts.txt"
printf 'Recovered FastOMA %s collection from %s\n' \
    "${FASTOMA_VERSION}" "${COLLECT_WORK_DIR}" > "${OUTDIR}/collection_recovery.txt"
