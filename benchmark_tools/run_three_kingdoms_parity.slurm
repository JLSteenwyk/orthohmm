#!/usr/bin/env bash
set -euo pipefail

: "${METHOD:?METHOD is required}"
: "${REPO_ROOT:?REPO_ROOT is required}"
: "${INPUT_DIR:?INPUT_DIR is required}"
: "${REFERENCE:?REFERENCE is required}"
: "${OUTDIR:?OUTDIR is required}"

readonly BASE=${REPO_ROOT}
readonly PYTHON_BIN=${PYTHON_BIN:-/home/bizon/anaconda3/bin/python}
readonly THREADS=${THREADS:-${SLURM_CPUS_PER_TASK:-32}}
readonly THREADS_PER_WORKER=${THREADS_PER_WORKER:-8}
readonly NORMALIZER="${BASE}/benchmark_tools/normalize_three_kingdoms_orthogroups.py"
readonly SCORER="${BASE}/three_kingdoms/score_against_busco.py"
readonly ALIGNER=${ALIGNER:-/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/SOFTWARE/mafft-7.525-with-extensions/bin/mafft}
readonly TREE_BUILDER=${TREE_BUILDER:-/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/SOFTWARE/FastTree_v220/FastTree}
readonly START_EPOCH=$(date +%s)

if [[ $(find "${INPUT_DIR}" -maxdepth 1 -type f -name '*.fasta' | wc -l) -ne 12 ]]; then
    echo "Expected 12 Three Kingdoms proteomes in ${INPUT_DIR}" >&2
    exit 1
fi
for required in "${REFERENCE}" "${NORMALIZER}" "${SCORER}"; do
    [[ -s "${required}" ]] || { echo "Missing required file: ${required}" >&2; exit 1; }
done
if [[ -e "${OUTDIR}/status.txt" || -e "${OUTDIR}/orthogroups.txt" ]]; then
    echo "Refusing to overwrite an existing parity result: ${OUTDIR}" >&2
    exit 1
fi

mkdir -p "${OUTDIR}"
printf 'running\n' > "${OUTDIR}/status.txt"
printf '%s\n' "${METHOD}" > "${OUTDIR}/method.txt"
git -C "${BASE}" rev-parse HEAD > "${OUTDIR}/source_commit.txt"
sha256sum "${INPUT_DIR}"/*.fasta > "${OUTDIR}/input.sha256"
sha256sum "${REFERENCE}" > "${OUTDIR}/reference.sha256"
printf 'slurm_job_id\t%s\ncpus\t%s\nstarted_epoch_s\t%s\n' \
    "${SLURM_JOB_ID:-local}" "${THREADS}" "${START_EPOCH}" \
    > "${OUTDIR}/run_metadata.tsv"

finish() {
    local exit_code=$?
    local end_epoch
    end_epoch=$(date +%s)
    printf '%s\n' "$((end_epoch - START_EPOCH))" > "${OUTDIR}/runner_wall_time_s.txt"
    printf 'finished_epoch_s\t%s\nexit_code\t%s\n' "${end_epoch}" "${exit_code}" \
        >> "${OUTDIR}/run_metadata.tsv"
    if [[ ${exit_code} -eq 0 ]]; then
        printf 'complete\n' > "${OUTDIR}/status.txt"
    else
        printf 'failed:%s\n' "${exit_code}" > "${OUTDIR}/status.txt"
    fi
}
trap finish EXIT

stage_fastas() {
    local stage_dir=$1
    mkdir -p "${stage_dir}"
    for source in "${INPUT_DIR}"/*.fasta; do
        cp -L --reflink=auto "${source}" "${stage_dir}/$(basename "${source}")"
    done
}

run_orthohmm_phylogeny() {
    local output_dir="${OUTDIR}/output"
    local metrics_json="${OUTDIR}/metrics.json"
    printf '%s\n' \
        'OrthoHMM high_sensitivity + inferred phylogeny satellite_v2' \
        > "${OUTDIR}/tool_version.txt"
    local inference_start
    inference_start=$(date +%s)
    "${PYTHON_BIN}" "${BASE}/benchmark_tools/benchmark_production.py" \
        "${INPUT_DIR}" "${output_dir}" "${metrics_json}" \
        --cpu "${THREADS}" \
        --threads-per-worker "${THREADS_PER_WORKER}" \
        --matrix BLOSUM62 \
        --evalue 1e-4 \
        --clustering leiden \
        --cpm-resolution 0.1 \
        --refinement-profile default \
        --accuracy-profile high_sensitivity \
        --phylogeny reconcile \
        --species-tree-mode infer \
        --species-tree-rooting min_variance \
        --phylogeny-candidates satellite_v2 \
        --phylogeny-root-rule species_overlap \
        --phylogeny-pair-rule positive_paralogy \
        --aligner "${ALIGNER}" \
        --tree-builder "${TREE_BUILDER}" \
        > "${OUTDIR}/run.log" 2>&1
    printf '%s\n' "$(( $(date +%s) - inference_start ))" \
        > "${OUTDIR}/inference_wall_time_s.txt"
    "${PYTHON_BIN}" "${NORMALIZER}" root-hogs \
        "${output_dir}/orthohmm_phylogeny/orthohmm_root_hogs.tsv" \
        "${OUTDIR}/orthogroups.txt" >> "${OUTDIR}/run.log"
    "${PYTHON_BIN}" - "${metrics_json}" "${OUTDIR}/peak_rss_kib.txt" <<'PY'
import json
import sys

metrics = json.load(open(sys.argv[1]))
with open(sys.argv[2], "w") as handle:
    handle.write(f"{int(metrics['peak_process_tree_rss_bytes']) // 1024}\n")
PY
}

run_orthofinder_full() {
    source "${BASE}/benchmark_tools/orthofinder.sh"
    local stage_dir="${OUTDIR}/run/input"
    stage_fastas "${stage_dir}"
    orthofinder_verify_install > "${OUTDIR}/tool_version.txt"
    local inference_start
    inference_start=$(date +%s)
    /usr/bin/time -v -o "${OUTDIR}/time.log" \
        "${ORTHOFINDER_BIN}" -f "${stage_dir}" -t "${THREADS}" -a "${THREADS}" -S diamond \
        > "${OUTDIR}/run.log" 2>&1
    printf '%s\n' "$(( $(date +%s) - inference_start ))" \
        > "${OUTDIR}/inference_wall_time_s.txt"
    local results_dir
    results_dir=$(orthofinder_results_dir "${stage_dir}")
    orthofinder_collect_orthogroups "${results_dir}" "${OUTDIR}"
    awk '/Maximum resident set size/ {print $NF}' "${OUTDIR}/time.log" \
        > "${OUTDIR}/peak_rss_kib.txt"
}

run_orthofinder_sequence_only() {
    : "${SOURCE_OUTDIR:?SOURCE_OUTDIR must point to the completed OrthoFinder 3.1.5 run}"
    source "${BASE}/benchmark_tools/orthofinder.sh"
    local stage_dir="${SOURCE_OUTDIR}/run/input"
    local results_dir
    results_dir=$(orthofinder_results_dir "${stage_dir}")
    local working_dir="${results_dir}/WorkingDirectory"
    local clusters="${working_dir}/clusters_OrthoFinder_I1.2.txt_id_pairs.txt"
    local sequence_ids="${working_dir}/SequenceIDs.txt"
    [[ -s "${clusters}" && -s "${sequence_ids}" ]] || {
        echo "Missing OrthoFinder sequence-only checkpoint" >&2
        exit 1
    }
    cp "${SOURCE_OUTDIR}/tool_version.txt" "${OUTDIR}/tool_version.txt"
    printf '%s\n' \
        'Derived from the matching full OrthoFinder 3.1.5 run at its MCL checkpoint.' \
        > "${OUTDIR}/runtime_note.txt"
    /usr/bin/time -v -o "${OUTDIR}/conversion_time.log" \
        "${PYTHON_BIN}" "${BASE}/benchmark_tools/orthofinder_mcl_to_orthogroups.py" \
        "${clusters}" "${sequence_ids}" "${OUTDIR}/orthogroups.txt" \
        > "${OUTDIR}/run.log" 2>&1
    "${PYTHON_BIN}" "${BASE}/benchmark_tools/derive_orthofinder_mcl_runtime.py" \
        "${OUTDIR}/inference_wall_time_s.txt" \
        "${results_dir}/Log.txt" "${SOURCE_OUTDIR}/run.log"
    sha256sum "${clusters}" "${sequence_ids}" > "${OUTDIR}/checkpoint.sha256"
}

run_fastoma() {
    : "${SOURCE_OUTDIR:?SOURCE_OUTDIR must point to the matching OrthoFinder run}"
    readonly FASTOMA_ROOT=${FASTOMA_ROOT:-/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/SOFTWARE/FastOMA}
    readonly FASTOMA_VERSION=${FASTOMA_VERSION:-0.3.5}
    readonly NEXTFLOW_BIN=${NEXTFLOW_BIN:-/home/bizon/bin/nextflow}
    readonly JAVA_HOME_17=${JAVA_HOME_17:-/home/bizon/anaconda3/pkgs/openjdk-17.0.18-ha668962_0/lib/jvm}
    readonly FASTOMA_CONFIG=${FASTOMA_CONFIG:-${BASE}/benchmark_tools/fastoma_qfo.config}
    readonly OMAMER_DB=${OMAMER_DB:-${BASE}/qfo_benchmark/data/LUCA.h5}
    readonly MAX_MEMORY=${MAX_MEMORY:-700.GB}
    local species_tree
    local species_trees=()
    mapfile -t species_trees < <(
        find "${SOURCE_OUTDIR}/run/input/OrthoFinder" -type f \
            -path '*/Species_Tree/SpeciesTree_rooted.txt' -print
    )
    if [[ ${#species_trees[@]} -ne 1 ]]; then
        echo "Expected one matching OrthoFinder species tree; found ${#species_trees[@]}" >&2
        exit 1
    fi
    species_tree=${species_trees[0]}
    local stage_dir="${OUTDIR}/run/input"
    local proteome_dir="${stage_dir}/proteome"
    local output_dir="${OUTDIR}/output"
    mkdir -p "${proteome_dir}"
    for source in "${INPUT_DIR}"/*.fasta; do
        cp -L --reflink=auto "${source}" \
            "${proteome_dir}/$(basename "${source}" .fasta).fa"
    done
    cp -L "${species_tree}" "${stage_dir}/species_tree.nwk"
    for required in "${NEXTFLOW_BIN}" "${FASTOMA_ROOT}/FastOMA.nf" \
        "${FASTOMA_CONFIG}" "${OMAMER_DB}"; do
        [[ -s "${required}" ]] || { echo "Missing FastOMA dependency: ${required}" >&2; exit 1; }
    done
    printf 'FastOMA version: %s\nFastOMA source: %s\n' \
        "${FASTOMA_VERSION}" "$(git -C "${FASTOMA_ROOT}" rev-parse HEAD)" \
        > "${OUTDIR}/tool_version.txt"
    sha256sum "${stage_dir}/species_tree.nwk" > "${OUTDIR}/species_tree.sha256"
    local inference_start
    inference_start=$(date +%s)
    cd "${FASTOMA_ROOT}"
    export JAVA_HOME="${JAVA_HOME_17}"
    export JAVA_CMD="${JAVA_HOME_17}/bin/java"
    /usr/bin/time -v -o "${OUTDIR}/time.log" \
        "${NEXTFLOW_BIN}" -c "${FASTOMA_CONFIG}" run FastOMA.nf \
        -profile docker -ansi-log false \
        --container_version "${FASTOMA_VERSION}" \
        --input_folder "${stage_dir}" \
        --output_folder "${output_dir}" \
        --omamer_db "${OMAMER_DB}" \
        --fasta_header_id_transformer noop \
        --force_pairwise_ortholog_generation true \
        --max_cpus "${THREADS}" \
        --max_memory "${MAX_MEMORY}" \
        -work-dir "${OUTDIR}/work" \
        > "${OUTDIR}/run.log" 2>&1
    printf '%s\n' "$(( $(date +%s) - inference_start ))" \
        > "${OUTDIR}/inference_wall_time_s.txt"
    [[ -s "${output_dir}/OrthologousGroups.tsv" ]] || {
        echo "FastOMA did not produce OrthologousGroups.tsv" >&2
        exit 1
    }
    "${PYTHON_BIN}" "${NORMALIZER}" fastoma \
        "${output_dir}/OrthologousGroups.tsv" "${OUTDIR}/orthogroups.txt" \
        >> "${OUTDIR}/run.log"
    awk '/Maximum resident set size/ {print $NF}' "${OUTDIR}/time.log" \
        > "${OUTDIR}/peak_rss_kib.txt"
}

run_orthomcl() {
    local workspace_root
    workspace_root=$(dirname "$(dirname "${BASE}")")
    readonly TEMPLATE_DIR=${ORTHOMCL_TEMPLATE_DIR:-${workspace_root}/SOFTWARE/ORTHOMCLV1.4}
    readonly ORTHOMCL_ENV=${ORTHOMCL_ENV:-${TEMPLATE_DIR}/venv_orthomcl}
    readonly ORTHOMCL_PERL=${ORTHOMCL_PERL:-${ORTHOMCL_ENV}/bin/perl}
    readonly TOOL_DIR="${OUTDIR}/orthomcl"
    readonly DATA_DIR="${TOOL_DIR}/data"
    [[ -x "${ORTHOMCL_PERL}" ]] || {
        echo "Missing OrthoMCL Perl interpreter: ${ORTHOMCL_PERL}" >&2
        exit 1
    }
    "${ORTHOMCL_PERL}" -MBio::SearchIO -e 1 || {
        echo "OrthoMCL Perl environment lacks Bio::SearchIO" >&2
        exit 1
    }
    mkdir -p "${DATA_DIR}"
    cp "${TEMPLATE_DIR}/orthomcl.pl" "${TEMPLATE_DIR}/orthomcl_module.pm" "${TOOL_DIR}/"
    "${PYTHON_BIN}" "${BASE}/benchmark_tools/configure_orthomcl_1_4.py" \
        "${TOOL_DIR}/orthomcl_module.pm" "${TOOL_DIR}" "${DATA_DIR}" "${THREADS}"
    local files=()
    for source in "${INPUT_DIR}"/*.fasta; do
        local target="$(basename "${source}" .fasta).fa"
        cp -L --reflink=auto "${source}" "${DATA_DIR}/${target}"
        files+=("${target}")
    done
    local file_list
    file_list=$(IFS=,; echo "${files[*]}")
    printf '%s\n' 'OrthoMCL 1.4 + MCL 02-063; 32-thread legacy BLAST' \
        > "${OUTDIR}/tool_version.txt"
    sha256sum "${TOOL_DIR}/orthomcl.pl" "${TOOL_DIR}/orthomcl_module.pm" \
        > "${OUTDIR}/tool.sha256"
    local inference_start
    inference_start=$(date +%s)
    cd "${TOOL_DIR}"
    /usr/bin/time -v -o "${OUTDIR}/time.log" \
        env PATH="${ORTHOMCL_ENV}/bin:${PATH}" PERL5LIB="${TOOL_DIR}:${PERL5LIB:-}" \
        "${ORTHOMCL_PERL}" "${TOOL_DIR}/orthomcl.pl" \
        --mode 1 --fa_files "${file_list}" \
        > "${OUTDIR}/run.log" 2>&1
    printf '%s\n' "$(( $(date +%s) - inference_start ))" \
        > "${OUTDIR}/inference_wall_time_s.txt"
    local native_output
    native_output=$(find "${TOOL_DIR}" -mindepth 2 -maxdepth 2 \
        -type f -name all_orthomcl.out -print -quit)
    [[ -s "${native_output}" ]] || { echo "Missing OrthoMCL output" >&2; exit 1; }
    cp "${native_output}" "${OUTDIR}/all_orthomcl.out"
    "${PYTHON_BIN}" "${NORMALIZER}" orthomcl \
        "${OUTDIR}/all_orthomcl.out" "${OUTDIR}/orthogroups.txt" \
        >> "${OUTDIR}/run.log"
    awk '/Maximum resident set size/ {print $NF}' "${OUTDIR}/time.log" \
        > "${OUTDIR}/peak_rss_kib.txt"
}

case "${METHOD}" in
    orthohmm_phylogeny_satellite_v2) run_orthohmm_phylogeny ;;
    orthofinder_3_1_5_full) run_orthofinder_full ;;
    orthofinder_3_1_5_sequence_only) run_orthofinder_sequence_only ;;
    fastoma_0_3_5) run_fastoma ;;
    orthomcl_1_4) run_orthomcl ;;
    *) echo "Unsupported Three Kingdoms method: ${METHOD}" >&2; exit 2 ;;
esac

"${PYTHON_BIN}" "${SCORER}" \
    --predictions "${OUTDIR}/orthogroups.txt" \
    --reference "${REFERENCE}" \
    --label "${METHOD}" \
    > "${OUTDIR}/score.txt"
sha256sum "${OUTDIR}/orthogroups.txt" "${OUTDIR}/score.txt" \
    > "${OUTDIR}/result.sha256"
