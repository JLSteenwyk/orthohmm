#!/usr/bin/env bash
set -euo pipefail

: "${REPO_ROOT:?REPO_ROOT is required}"
: "${INPUT_DIR:?INPUT_DIR is required}"
: "${OUTDIR:?OUTDIR is required}"

readonly BASE=${REPO_ROOT}
readonly SCRIPT_DIR="${BASE}/benchmark_tools"
readonly PYTHON_BIN=${PYTHON_BIN:-/home/bizon/anaconda3/bin/python}
readonly THREADS=${THREADS:-${SLURM_CPUS_PER_TASK:-180}}
readonly PAIR_WORKERS=${ORTHOMCL_PAIR_WORKERS:-64}
readonly EXPECTED_TAXA=${EXPECTED_TAXA:-78}
readonly EXPECTED_SEQUENCES=${EXPECTED_SEQUENCES:-976504}
readonly TEMPLATE_DIR=${ORTHOMCL_TEMPLATE_DIR:-/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/SOFTWARE/ORTHOMCLV1.4}
readonly ORTHOMCL_ENV=${ORTHOMCL_ENV:-${TEMPLATE_DIR}/venv_orthomcl}
readonly ORTHOMCL_PERL=${ORTHOMCL_PERL:-${ORTHOMCL_ENV}/bin/perl}
readonly BLASTALL=${BLASTALL:-/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/SOFTWARE/blast-2.2.13/bin/blastall}
readonly FORMATDB=${FORMATDB:-/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/SOFTWARE/blast-2.2.13/bin/formatdb}
readonly MCL=${MCL:-/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/SOFTWARE/mcl-02-063/bin/mcl}
readonly TOOL_DIR="${OUTDIR}/orthomcl"
readonly WORK_DIR="${OUTDIR}/work"
readonly STAGE_DIR="${OUTDIR}/stages"
readonly ALL_FASTA="${WORK_DIR}/all.fa"
readonly GG_FILE="${WORK_DIR}/all.gg"
readonly BLAST_FILE="${WORK_DIR}/all.blast"
readonly BPO_FILE="${WORK_DIR}/all.bpo"
readonly START_EPOCH=$(date +%s)

for required in "${PYTHON_BIN}" "${ORTHOMCL_PERL}" "${BLASTALL}" "${FORMATDB}" \
    "${MCL}" "${TEMPLATE_DIR}/orthomcl.pl" "${TEMPLATE_DIR}/orthomcl_module.pm" \
    "${SCRIPT_DIR}/prepare_orthomcl_inputs.py" \
    "${SCRIPT_DIR}/configure_orthomcl_1_4.py" \
    "${SCRIPT_DIR}/convert_orthomcl_blast.py" \
    "${SCRIPT_DIR}/parallelize_orthomcl_pairs.py" \
    "${SCRIPT_DIR}/normalize_three_kingdoms_orthogroups.py" \
    "${SCRIPT_DIR}/orthomcl_matrix_to_pairwise.py" \
    "${SCRIPT_DIR}/qfo_filter_pairs.py" \
    "${BASE}/qfo_benchmark/benchmark-webservice/reference_data/2020/mapping.json.gz"; do
    [[ -s "${required}" && -x "${required}" || -s "${required}" ]] || {
        echo "Missing required file: ${required}" >&2
        exit 1
    }
done
[[ ${THREADS} =~ ^[1-9][0-9]*$ ]] || { echo "Invalid THREADS: ${THREADS}" >&2; exit 2; }
[[ ${PAIR_WORKERS} =~ ^[1-9][0-9]*$ ]] || {
    echo "Invalid ORTHOMCL_PAIR_WORKERS: ${PAIR_WORKERS}" >&2
    exit 2
}
[[ $(find "${INPUT_DIR}" -maxdepth 1 -type f -name '*.fasta' | wc -l) -eq ${EXPECTED_TAXA} ]] || {
    echo "Expected ${EXPECTED_TAXA} QfO proteomes in ${INPUT_DIR}" >&2
    exit 1
}

mkdir -p "${OUTDIR}" "${WORK_DIR}" "${STAGE_DIR}"
if [[ -s "${OUTDIR}/status.txt" && $(<"${OUTDIR}/status.txt") == complete ]]; then
    echo "OrthoMCL QfO result is already complete: ${OUTDIR}"
    exit 0
fi

current_manifest=$(mktemp "${OUTDIR}/input.sha256.XXXXXX")
sha256sum "${INPUT_DIR}"/*.fasta > "${current_manifest}"
if [[ -s "${OUTDIR}/input.sha256" ]]; then
    cmp --silent "${OUTDIR}/input.sha256" "${current_manifest}" || {
        echo "QfO input manifest changed; refusing to combine incompatible stages" >&2
        rm -f "${current_manifest}"
        exit 1
    }
    rm -f "${current_manifest}"
else
    mv "${current_manifest}" "${OUTDIR}/input.sha256"
fi

printf 'running:initialization\n' > "${OUTDIR}/status.txt"
printf '%s\n' orthomcl_1_4 > "${OUTDIR}/method.txt"
git -C "${BASE}" rev-parse HEAD > "${OUTDIR}/source_commit.txt"
printf '%s\n' \
    "OrthoMCL 1.4 + legacy NCBI BLAST 2.2.13 (${THREADS} threads) + MCL 02-063; exact-compatible streaming BPO conversion; ${PAIR_WORKERS} pair workers" \
    > "${OUTDIR}/tool_version.txt"
printf 'run_started_epoch_s\t%s\nslurm_job_id\t%s\ncpus\t%s\npair_workers\t%s\n' \
    "${START_EPOCH}" "${SLURM_JOB_ID:-local}" "${THREADS}" "${PAIR_WORKERS}" \
    >> "${OUTDIR}/run_metadata.tsv"

finish() {
    local exit_code=$?
    local end_epoch
    end_epoch=$(date +%s)
    printf 'run_finished_epoch_s\t%s\nexit_code\t%s\nrunner_wall_time_s\t%s\n' \
        "${end_epoch}" "${exit_code}" "$((end_epoch - START_EPOCH))" \
        >> "${OUTDIR}/run_metadata.tsv"
    if [[ ${exit_code} -eq 0 && -s "${OUTDIR}/pairs.qfo.tsv" ]]; then
        printf 'complete\n' > "${OUTDIR}/status.txt"
    else
        printf 'failed:%s\n' "${exit_code}" > "${OUTDIR}/status.txt"
    fi
}
trap finish EXIT

stage_begin() {
    CURRENT_STAGE=$1
    STAGE_START_EPOCH=$(date +%s)
    printf 'running:%s\n' "${CURRENT_STAGE}" > "${OUTDIR}/status.txt"
    printf '%s\n' "${STAGE_START_EPOCH}" > "${STAGE_DIR}/${CURRENT_STAGE}.started_epoch_s"
}

stage_complete() {
    local end_epoch
    end_epoch=$(date +%s)
    printf '%s\t%s\n' "${CURRENT_STAGE}" "$((end_epoch - STAGE_START_EPOCH))" \
        >> "${OUTDIR}/stage_wall_times.tsv"
    touch "${STAGE_DIR}/${CURRENT_STAGE}.complete"
}

if [[ ! -e "${STAGE_DIR}/tools.complete" ]]; then
    stage_begin tools
    rm -rf "${TOOL_DIR}"
    mkdir -p "${TOOL_DIR}/data"
    cp "${TEMPLATE_DIR}/orthomcl.pl" "${TEMPLATE_DIR}/orthomcl_module.pm" "${TOOL_DIR}/"
    "${PYTHON_BIN}" "${SCRIPT_DIR}/configure_orthomcl_1_4.py" \
        "${TOOL_DIR}/orthomcl_module.pm" "${TOOL_DIR}" "${TOOL_DIR}/data" "${THREADS}"
    "${PYTHON_BIN}" "${SCRIPT_DIR}/parallelize_orthomcl_pairs.py" \
        "${TOOL_DIR}/orthomcl.pl"
    "${ORTHOMCL_PERL}" -I "${TOOL_DIR}" -c "${TOOL_DIR}/orthomcl.pl"
    sha256sum "${TOOL_DIR}/orthomcl.pl" "${TOOL_DIR}/orthomcl_module.pm" \
        > "${OUTDIR}/tool.sha256"
    stage_complete
fi
sha256sum --check --status "${OUTDIR}/tool.sha256" || {
    echo "Configured OrthoMCL tools do not match their recorded checksums" >&2
    exit 1
}

if [[ ! -e "${STAGE_DIR}/prepare.complete" ]]; then
    stage_begin prepare
    "${PYTHON_BIN}" "${SCRIPT_DIR}/prepare_orthomcl_inputs.py" \
        "${INPUT_DIR}" "${ALL_FASTA}" "${GG_FILE}" \
        --summary "${OUTDIR}/input_summary.json" \
        > "${OUTDIR}/prepare.log" 2>&1
    "${PYTHON_BIN}" - "${OUTDIR}/input_summary.json" \
        "${EXPECTED_TAXA}" "${EXPECTED_SEQUENCES}" <<'PY'
import json
import sys

summary = json.load(open(sys.argv[1]))
expected = (int(sys.argv[2]), int(sys.argv[3]))
observed = (summary["taxon_count"], summary["sequence_count"])
if observed != expected:
    raise SystemExit(f"Unexpected QfO input dimensions: {observed} != {expected}")
PY
    stage_complete
fi
[[ -s "${ALL_FASTA}" && -s "${GG_FILE}" ]] || {
    echo "Prepared OrthoMCL inputs are missing" >&2
    exit 1
}

if [[ ! -e "${STAGE_DIR}/formatdb.complete" ]]; then
    stage_begin formatdb
    rm -f "${ALL_FASTA}.p"{hr,in,sq}
    /usr/bin/time -v -o "${OUTDIR}/formatdb.time.log" \
        "${FORMATDB}" -i "${ALL_FASTA}" -p t \
        > "${OUTDIR}/formatdb.log" 2>&1
    for database_file in "${ALL_FASTA}.p"{hr,in,sq}; do
        [[ -s "${database_file}" ]] || { echo "Missing BLAST database file: ${database_file}" >&2; exit 1; }
    done
    stage_complete
fi
for database_file in "${ALL_FASTA}.p"{hr,in,sq}; do
    [[ -s "${database_file}" ]] || { echo "Missing BLAST database file: ${database_file}" >&2; exit 1; }
done

if [[ ! -e "${STAGE_DIR}/blast.complete" ]]; then
    stage_begin blast
    rm -f "${BLAST_FILE}.partial"
    /usr/bin/time -v -o "${OUTDIR}/blast.time.log" \
        "${BLASTALL}" -p blastp -i "${ALL_FASTA}" -d "${ALL_FASTA}" \
        -e 1e-5 -o "${BLAST_FILE}.partial" -m 8 -a "${THREADS}" -v 1000 -b 1000 \
        > "${OUTDIR}/blast.log" 2>&1
    [[ -s "${BLAST_FILE}.partial" ]] || { echo "BLAST produced no hits" >&2; exit 1; }
    mv "${BLAST_FILE}.partial" "${BLAST_FILE}"
    stage_complete
fi
[[ -s "${BLAST_FILE}" ]] || { echo "BLAST result is missing" >&2; exit 1; }

if [[ ! -e "${STAGE_DIR}/bpo.complete" ]]; then
    stage_begin bpo
    rm -f "${BPO_FILE}" "${BPO_FILE}.partial"
    /usr/bin/time -v -o "${OUTDIR}/bpo.time.log" \
        "${PYTHON_BIN}" "${SCRIPT_DIR}/convert_orthomcl_blast.py" \
        "${BLAST_FILE}" "${ALL_FASTA}" "${BPO_FILE}.partial" \
        > "${OUTDIR}/bpo.log" 2>&1
    [[ -s "${BPO_FILE}.partial" ]] || { echo "BPO conversion produced no hits" >&2; exit 1; }
    mv "${BPO_FILE}.partial" "${BPO_FILE}"
    stage_complete
fi
[[ -s "${BPO_FILE}" ]] || { echo "BPO result is missing" >&2; exit 1; }

readonly BPO_INDEX="${WORK_DIR}/all_bpo.idx"
readonly BPO_SE="${WORK_DIR}/all_bpo.se"
if [[ ! -e "${STAGE_DIR}/indexes.complete" ]]; then
    stage_begin indexes
    rm -f "${BPO_INDEX}" "${BPO_SE}" "${BPO_INDEX}.partial" "${BPO_SE}.partial"
    /usr/bin/time -v -o "${OUTDIR}/indexes.time.log" \
        env PERL5LIB="${TOOL_DIR}:${PERL5LIB:-}" \
        "${ORTHOMCL_PERL}" -I "${TOOL_DIR}" -Morthomcl_module -e \
        'constructIDX_for_bpofile($ARGV[0], $ARGV[1]); constructSE_for_bpofile($ARGV[0], $ARGV[2]);' \
        "${BPO_FILE}" "${BPO_INDEX}.partial" "${BPO_SE}.partial" \
        > "${OUTDIR}/indexes.log" 2>&1
    [[ -s "${BPO_INDEX}.partial" && -s "${BPO_SE}.partial" ]] || {
        echo "OrthoMCL BPO indexes are missing" >&2
        exit 1
    }
    mv "${BPO_INDEX}.partial" "${BPO_INDEX}"
    mv "${BPO_SE}.partial" "${BPO_SE}"
    stage_complete
fi
[[ -s "${BPO_INDEX}" && -s "${BPO_SE}" ]] || {
    echo "OrthoMCL BPO indexes are missing" >&2
    exit 1
}

if [[ ! -e "${STAGE_DIR}/inference.complete" ]]; then
    stage_begin inference
    mapfile -t existing_outputs < <(
        find "${TOOL_DIR}" -mindepth 2 -maxdepth 2 -type f \
            -name all_orthomcl.out -size +0c -print
    )
    if [[ ${#existing_outputs[@]} -eq 0 ]]; then
        cd "${TOOL_DIR}"
        /usr/bin/time -v -o "${OUTDIR}/inference.time.log" \
            env PATH="${ORTHOMCL_ENV}/bin:${PATH}" \
            PERL5LIB="${TOOL_DIR}:${PERL5LIB:-}" \
            ORTHOMCL_PAIR_WORKERS="${PAIR_WORKERS}" \
            "${ORTHOMCL_PERL}" "${TOOL_DIR}/orthomcl.pl" \
            --mode 4 --bpo_file "${BPO_FILE}" --gg_file "${GG_FILE}" \
            > "${OUTDIR}/inference.log" 2>&1
        cd "${BASE}"
        mapfile -t existing_outputs < <(
            find "${TOOL_DIR}" -mindepth 2 -maxdepth 2 -type f \
                -name all_orthomcl.out -size +0c -print
        )
    fi
    if [[ ${#existing_outputs[@]} -ne 1 ]]; then
        echo "Expected one completed native OrthoMCL result; found ${#existing_outputs[@]}" >&2
        exit 1
    fi
    cp "${existing_outputs[0]}" "${OUTDIR}/all_orthomcl.out.partial"
    mv "${OUTDIR}/all_orthomcl.out.partial" "${OUTDIR}/all_orthomcl.out"
    stage_complete
fi
[[ -s "${OUTDIR}/all_orthomcl.out" ]] || {
    echo "Native OrthoMCL result is missing" >&2
    exit 1
}
awk '$1 ~ /^(formatdb|blast|bpo|indexes|inference)$/ {total += $2} END {print total}' \
    "${OUTDIR}/stage_wall_times.tsv" > "${OUTDIR}/inference_wall_time_s.txt"

if [[ ! -e "${STAGE_DIR}/pairs.complete" ]]; then
    stage_begin pairs
    rm -f "${OUTDIR}/orthogroups.txt.partial" "${OUTDIR}/pairs.tsv.partial" \
        "${OUTDIR}/pairs.qfo.tsv.partial"
    "${PYTHON_BIN}" "${SCRIPT_DIR}/normalize_three_kingdoms_orthogroups.py" \
        orthomcl "${OUTDIR}/all_orthomcl.out" "${OUTDIR}/orthogroups.txt.partial" \
        > "${OUTDIR}/pairs.log" 2>&1
    mv "${OUTDIR}/orthogroups.txt.partial" "${OUTDIR}/orthogroups.txt"
    mapfile -t native_outputs < <(
        find "${TOOL_DIR}" -mindepth 2 -maxdepth 2 -type f \
            -name all_orthomcl.out -size +0c -print
    )
    if [[ ${#native_outputs[@]} -ne 1 ]]; then
        echo "Expected one completed native OrthoMCL result; found ${#native_outputs[@]}" >&2
        exit 1
    fi
    readonly NATIVE_RUN_DIR=$(dirname "${native_outputs[0]}")
    readonly NATIVE_MATRIX="${NATIVE_RUN_DIR}/tmp/all_ortho.mtx"
    readonly NATIVE_INDEX="${NATIVE_RUN_DIR}/tmp/all_ortho.idx"
    [[ -s "${NATIVE_MATRIX}" && -s "${NATIVE_INDEX}" ]] || {
        echo "Native OrthoMCL matrix or index is missing from ${NATIVE_RUN_DIR}" >&2
        exit 1
    }
    "${PYTHON_BIN}" "${SCRIPT_DIR}/orthomcl_matrix_to_pairwise.py" \
        "${NATIVE_MATRIX}" "${NATIVE_INDEX}" "${GG_FILE}" \
        > "${OUTDIR}/pairs.tsv.partial" 2>> "${OUTDIR}/pairs.log"
    mv "${OUTDIR}/pairs.tsv.partial" "${OUTDIR}/pairs.tsv"
    "${PYTHON_BIN}" "${SCRIPT_DIR}/qfo_filter_pairs.py" \
        "${OUTDIR}/pairs.tsv" \
        "${BASE}/qfo_benchmark/benchmark-webservice/reference_data/2020/mapping.json.gz" \
        "${OUTDIR}/pairs.qfo.tsv.partial" >> "${OUTDIR}/pairs.log" 2>&1
    mv "${OUTDIR}/pairs.qfo.tsv.partial" "${OUTDIR}/pairs.qfo.tsv"
    [[ -s "${OUTDIR}/pairs.qfo.tsv" ]] || { echo "No QfO pairs were retained" >&2; exit 1; }
    printf '%s\n' "native OrthoMCL cross-species ortholog/co-ortholog matrix edges" \
        > "${OUTDIR}/pair_source.txt"
    sha256sum "${OUTDIR}/all_orthomcl.out" "${OUTDIR}/orthogroups.txt" \
        "${NATIVE_MATRIX}" "${NATIVE_INDEX}" "${OUTDIR}/pairs.qfo.tsv" \
        > "${OUTDIR}/result.sha256"
    stage_complete
fi
[[ -s "${OUTDIR}/orthogroups.txt" && -s "${OUTDIR}/pairs.qfo.tsv" ]] || {
    echo "OrthoMCL QfO benchmark artifacts are incomplete" >&2
    exit 1
}
