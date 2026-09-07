#!/usr/bin/env bash
set -euo pipefail

readonly SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
readonly BASE=$(dirname "${SCRIPT_DIR}")
readonly INPUT_DIR=${INPUT_DIR:-${BASE}/three_kingdoms/input}
readonly REFERENCE=${REFERENCE:-${BASE}/three_kingdoms/busco/reference_orthogroups.txt}
readonly RUN_ID=${RUN_ID:-parity_20260907}
readonly RUN_ROOT=${RUN_ROOT:-${BASE}/three_kingdoms/results/${RUN_ID}}
readonly LOG_DIR="${RUN_ROOT}/slurm_logs"
readonly RUNNER="${SCRIPT_DIR}/run_three_kingdoms_parity.slurm"

mkdir -p "${LOG_DIR}"

submit() {
    local method=$1
    local cpus=$2
    local mem=$3
    local dependency=${4:-}
    local extra_export=${5:-}
    local outdir="${RUN_ROOT}/${method}"
    if [[ -s "${outdir}/status.txt" ]] && grep -Fxq complete "${outdir}/status.txt"; then
        printf 'complete:%s\n' "${method}"
        return
    fi
    if [[ -e "${outdir}/status.txt" ]]; then
        echo "Existing incomplete result requires review: ${outdir}" >&2
        exit 1
    fi
    local dependency_arg=()
    [[ -n "${dependency}" ]] && dependency_arg=(--dependency="${dependency}")
    local export_value="ALL,METHOD=${method},REPO_ROOT=${BASE},INPUT_DIR=${INPUT_DIR},REFERENCE=${REFERENCE},OUTDIR=${outdir}"
    [[ -n "${extra_export}" ]] && export_value+=",${extra_export}"
    sbatch --parsable \
        --job-name="tk_${method}" \
        --time=14-00:00:00 \
        --cpus-per-task="${cpus}" \
        --mem="${mem}" \
        --output="${LOG_DIR}/${method}_%j.out" \
        --error="${LOG_DIR}/${method}_%j.err" \
        "${dependency_arg[@]}" \
        --export="${export_value}" \
        "${RUNNER}"
}

of_job=$(submit orthofinder_3_1_5_full 32 250G)
phylogeny_job=$(submit orthohmm_phylogeny_satellite_v2 32 250G)

if [[ ${of_job} == complete:* ]]; then
    echo "This submitter expects a new OrthoFinder run when launching dependencies." >&2
    exit 1
fi
sequence_job=$(submit orthofinder_3_1_5_sequence_only 1 32G \
    "afterok:${of_job}" \
    "SOURCE_OUTDIR=${RUN_ROOT}/orthofinder_3_1_5_full")
fastoma_job=$(submit fastoma_0_3_5 128 700G \
    "afterok:${of_job}" \
    "SOURCE_OUTDIR=${RUN_ROOT}/orthofinder_3_1_5_full")

orthomcl_dependency="afterany:${phylogeny_job}:${fastoma_job}"
orthomcl_job=$(submit orthomcl_1_4 32 700G "${orthomcl_dependency}")

printf 'OrthoFinder 3.1.5 full:          %s\n' "${of_job}"
printf 'OrthoFinder 3.1.5 sequence-only: %s\n' "${sequence_job}"
printf 'OrthoHMM phylogeny satellite_v2: %s\n' "${phylogeny_job}"
printf 'FastOMA 0.3.5:                   %s\n' "${fastoma_job}"
printf 'OrthoMCL 1.4:                    %s\n' "${orthomcl_job}"
