#!/usr/bin/env bash
set -euo pipefail

: "${REPO_ROOT:?REPO_ROOT is required}"
: "${REFERENCE:?REFERENCE is required}"
: "${OUTDIR:?OUTDIR is required}"
: "${BPO_FILE:?BPO_FILE is required}"
: "${BLAST_FILE:?BLAST_FILE is required}"
: "${GG_FILE:?GG_FILE is required}"
: "${BPO_CONVERTER_WALL_S:?BPO_CONVERTER_WALL_S is required}"

readonly BASE=${REPO_ROOT}
readonly PYTHON_BIN=${PYTHON_BIN:-/home/bizon/anaconda3/bin/python}
readonly THREADS=${THREADS:-${SLURM_CPUS_PER_TASK:-32}}
readonly TOOL_DIR="${OUTDIR}/orthomcl"
readonly ORTHOMCL_ENV=${ORTHOMCL_ENV:-/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/SOFTWARE/ORTHOMCLV1.4/venv_orthomcl}
readonly ORTHOMCL_PERL=${ORTHOMCL_PERL:-${ORTHOMCL_ENV}/bin/perl}
readonly NORMALIZER="${BASE}/benchmark_tools/normalize_three_kingdoms_orthogroups.py"
readonly SCORER="${BASE}/three_kingdoms/score_against_busco.py"
readonly PAIR_PARALLELIZER="${BASE}/benchmark_tools/parallelize_orthomcl_pairs.py"
readonly PAIR_WORKERS=${ORTHOMCL_PAIR_WORKERS:-1}
readonly INDEX_WALL_S=${ORTHOMCL_INDEX_WALL_S:-0}
readonly RECOVERY_START_EPOCH=$(date +%s)

for required in "${BPO_FILE}" "${BLAST_FILE}" "${GG_FILE}" "${REFERENCE}" \
    "${TOOL_DIR}/orthomcl.pl" "${TOOL_DIR}/orthomcl_module.pm" \
    "${NORMALIZER}" "${SCORER}" "${PAIR_PARALLELIZER}"; do
    [[ -s "${required}" ]] || { echo "Missing recovery input: ${required}" >&2; exit 1; }
done
[[ -x "${ORTHOMCL_PERL}" ]] || { echo "Missing Perl: ${ORTHOMCL_PERL}" >&2; exit 1; }
"${ORTHOMCL_PERL}" -MBio::SearchIO -e 1
[[ -s "${OUTDIR}/bpo_prefix_validation.txt" ]] || {
    echo "Missing native-prefix validation for ${BPO_FILE}" >&2
    exit 1
}
[[ ! -e "${OUTDIR}/orthogroups.txt" ]] || {
    echo "Refusing to overwrite an existing OrthoMCL result" >&2
    exit 1
}

original_start_epoch=$(awk -F '\t' '$1 == "started_epoch_s" {print $2; exit}' \
    "${OUTDIR}/run_metadata.tsv")
[[ ${original_start_epoch} =~ ^[0-9]+$ ]] || {
    echo "Invalid original start epoch: ${original_start_epoch}" >&2
    exit 1
}
blast_end_epoch=$(stat -c %Y "${BLAST_FILE}")
blast_wall_s=$((blast_end_epoch - original_start_epoch))
[[ ${blast_wall_s} -gt 0 ]] || { echo "Invalid BLAST stage duration" >&2; exit 1; }

finish() {
    local exit_code=$?
    local end_epoch
    end_epoch=$(date +%s)
    printf '%s\n' "$((end_epoch - RECOVERY_START_EPOCH))" \
        > "${OUTDIR}/recovery_runner_wall_time_s.txt"
    printf 'recovery_finished_epoch_s\t%s\nrecovery_exit_code\t%s\n' \
        "${end_epoch}" "${exit_code}" >> "${OUTDIR}/run_metadata.tsv"
    if [[ ${exit_code} -eq 0 && -s "${OUTDIR}/orthogroups.txt" && \
          -s "${OUTDIR}/score.txt" ]]; then
        printf 'complete\n' > "${OUTDIR}/status.txt"
    else
        printf 'failed_recovery:%s:incomplete\n' "${exit_code}" \
            > "${OUTDIR}/status.txt"
    fi
}
trap finish EXIT

printf 'running_recovery\n' > "${OUTDIR}/status.txt"
printf 'recovery_slurm_job_id\t%s\nrecovery_started_epoch_s\t%s\n' \
    "${SLURM_JOB_ID:-local}" "${RECOVERY_START_EPOCH}" \
    >> "${OUTDIR}/run_metadata.tsv"
git -C "${BASE}" rev-parse HEAD > "${OUTDIR}/source_commit.txt"
printf '%s\n' \
    "OrthoMCL 1.4 + MCL 02-063; 32-thread legacy BLAST; exact-compatible BPO conversion; ${PAIR_WORKERS} pair workers" \
    > "${OUTDIR}/tool_version.txt"
if [[ ${PAIR_WORKERS} -gt 1 ]]; then
    if ! grep -q 'ORTHOMCL_PAIR_WORKERS' "${TOOL_DIR}/orthomcl.pl"; then
        "${PYTHON_BIN}" "${PAIR_PARALLELIZER}" "${TOOL_DIR}/orthomcl.pl"
    fi
    "${ORTHOMCL_PERL}" -I "${TOOL_DIR}" -c "${TOOL_DIR}/orthomcl.pl"
fi
sha256sum "${TOOL_DIR}/orthomcl.pl" "${TOOL_DIR}/orthomcl_module.pm" \
    > "${OUTDIR}/recovery_tool.sha256"
sha256sum "${BLAST_FILE}" "${BPO_FILE}" "${GG_FILE}" \
    > "${OUTDIR}/recovery_inputs.sha256"
printf '%s\n' "${blast_wall_s}" > "${OUTDIR}/blast_wall_time_s.txt"
printf '%s\n' "${BPO_CONVERTER_WALL_S}" \
    > "${OUTDIR}/bpo_conversion_wall_time_s.txt"
printf '%s\n' "${INDEX_WALL_S}" > "${OUTDIR}/bpo_index_wall_time_s.txt"

bpo_base=${BPO_FILE%.*}
if [[ ${ORTHOMCL_REUSE_INDEXES:-0} == 1 ]]; then
    for index in "${bpo_base}_bpo.idx" "${bpo_base}_bpo.se"; do
        [[ -s "${index}" ]] || { echo "Missing reusable index: ${index}" >&2; exit 1; }
    done
else
    for stale_index in "${bpo_base}_bpo.idx" "${bpo_base}_bpo.se"; do
        [[ ! -e "${stale_index}" ]] || {
            echo "Refusing to reuse an unverified OrthoMCL index: ${stale_index}" >&2
            exit 1
        }
    done
fi

downstream_start_epoch=$(date +%s)
cd "${TOOL_DIR}"
/usr/bin/time -v -o "${OUTDIR}/downstream_time.log" \
    env PATH="${ORTHOMCL_ENV}/bin:${PATH}" PERL5LIB="${TOOL_DIR}:${PERL5LIB:-}" \
    ORTHOMCL_PAIR_WORKERS="${PAIR_WORKERS}" \
    "${ORTHOMCL_PERL}" "${TOOL_DIR}/orthomcl.pl" \
    --mode 4 --bpo_file "${BPO_FILE}" --gg_file "${GG_FILE}" \
    > "${OUTDIR}/recovery_run.log" 2>&1
downstream_wall_s=$(( $(date +%s) - downstream_start_epoch ))
printf '%s\n' "${downstream_wall_s}" > "${OUTDIR}/downstream_wall_time_s.txt"
printf '%s\n' \
    "$((blast_wall_s + BPO_CONVERTER_WALL_S + INDEX_WALL_S + downstream_wall_s))" \
    > "${OUTDIR}/inference_wall_time_s.txt"

mapfile -t native_outputs < <(
    find "${TOOL_DIR}" -mindepth 2 -maxdepth 2 -type f \
        -name all_orthomcl.out -print
)
if [[ ${#native_outputs[@]} -ne 1 ]]; then
    echo "Expected one completed OrthoMCL output; found ${#native_outputs[@]}" >&2
    exit 1
fi
[[ -s ${native_outputs[0]} ]] || { echo "Empty OrthoMCL output" >&2; exit 1; }
cp "${native_outputs[0]}" "${OUTDIR}/all_orthomcl.out"
"${PYTHON_BIN}" "${NORMALIZER}" orthomcl \
    "${OUTDIR}/all_orthomcl.out" "${OUTDIR}/orthogroups.txt" \
    >> "${OUTDIR}/recovery_run.log"
"${PYTHON_BIN}" "${SCORER}" \
    --predictions "${OUTDIR}/orthogroups.txt" \
    --reference "${REFERENCE}" \
    --label orthomcl_1_4 > "${OUTDIR}/score.txt"

peak_rss_kib=0
for time_log in "${OUTDIR}/time.log" "${OUTDIR}/bpo_converter.time.log" \
    "${OUTDIR}/downstream_time.log"; do
    [[ -s "${time_log}" ]] || continue
    value=$(awk '/Maximum resident set size/ {print $NF}' "${time_log}")
    if [[ ${value:-0} =~ ^[0-9]+$ && ${value} -gt ${peak_rss_kib} ]]; then
        peak_rss_kib=${value}
    fi
done
printf '%s\n' "${peak_rss_kib}" > "${OUTDIR}/peak_rss_kib.txt"
sha256sum "${OUTDIR}/all_orthomcl.out" "${OUTDIR}/orthogroups.txt" \
    "${OUTDIR}/score.txt" > "${OUTDIR}/result.sha256"
