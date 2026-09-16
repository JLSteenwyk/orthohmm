#!/usr/bin/env bash
set -euo pipefail

: "${REPO_ROOT:?REPO_ROOT is required}"
readonly BASE=${REPO_ROOT}
readonly INPUT_DIR=${INPUT_DIR:-${BASE}/qfo_benchmark/input}
readonly ORIGINAL=${BASE}/qfo_benchmark/results/orthomcl_1_4
readonly OUT=${BASE}/qfo_benchmark/results/orthomcl_1_4_final_groups
readonly SCORING=${BASE}/qfo_benchmark/scoring

# New namespace only: never overwrite the completed pre-clustering assessment.
for path in "${OUT}" "${SCORING}/orthomcl_1_4_final_groups" \
    "${SCORING}/work-omcl_final" "${SCORING}/launch-orthomcl_1_4_final_groups"; do
    if [[ -e "${path}" ]]; then
        echo "Refusing to overwrite existing artifact: ${path}" >&2
        exit 2
    fi
done
sha256sum --check --status "${ORIGINAL}/input.sha256"
sha256sum --check --status "${ORIGINAL}/result.sha256"
mkdir "${OUT}"
trap 'printf "exit_code\t%s\n" "$?" >> "${OUT}/run_metadata.tsv"' EXIT
git -C "${BASE}" rev-parse HEAD > "${OUT}/conversion_source_commit.txt"
printf 'job_id\t%s\nstarted\t%s\n' "${SLURM_JOB_ID:-local}" "$(date -u +%FT%TZ)" \
    > "${OUT}/run_metadata.tsv"
printf '%s\n' 'Cross-species clique expansion of final MCL groups; no new inference or parameter selection.' \
    > "${OUT}/pair_source.txt"
sha256sum "${ORIGINAL}/all_orthomcl.out" "${ORIGINAL}/orthogroups.txt" \
    "${BASE}/qfo_benchmark/og_to_pairwise.py" > "${OUT}/inputs.sha256"
/usr/bin/time -v -o "${OUT}/conversion.time.log" \
    python "${BASE}/qfo_benchmark/og_to_pairwise.py" \
    "${ORIGINAL}/orthogroups.txt" "${INPUT_DIR}" \
    > "${OUT}/pairs.tsv.partial" 2> "${OUT}/conversion.log"
mv "${OUT}/pairs.tsv.partial" "${OUT}/pairs.tsv"
export METHOD=orthomcl_1_4_final_groups PAIRS="${OUT}/pairs.tsv"
/usr/bin/time -v -o "${OUT}/scoring.time.log" \
    bash "${BASE}/benchmark_tools/run_qfo_scoring.slurm"
sha256sum "${OUT}/pairs.tsv" "${OUT}/pairs.qfo.tsv" > "${OUT}/result.sha256"
printf 'finished\t%s\n' "$(date -u +%FT%TZ)" >> "${OUT}/run_metadata.tsv"
