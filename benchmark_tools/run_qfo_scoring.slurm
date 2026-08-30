#!/usr/bin/env bash
set -euo pipefail

: "${METHOD:?METHOD must name the QfO result}"
: "${PAIRS:?PAIRS must point to the QfO pair predictions}"

case "${METHOD}" in
    orthofinder_v3_diamond|orthofinder_v3_mmseqs|orthofinder_v3_sequence_only|orthohmm_high_sensitivity|orthohmm_phylogeny_inferred) ;;
    *)
        echo "Unsupported scoring method: ${METHOD}" >&2
        exit 2
        ;;
esac

if [[ ! -s "${PAIRS}" ]]; then
    echo "Pair predictions are missing or empty: ${PAIRS}" >&2
    exit 1
fi

readonly BASE=${REPO_ROOT:-$(dirname "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)")}
readonly SCRIPT_DIR="${BASE}/benchmark_tools"
readonly PIPELINE=${QFO_PIPELINE:-${BASE}/qfo_benchmark/benchmark-webservice}
readonly SCORING_ROOT=${QFO_SCORING_ROOT:-${BASE}/qfo_benchmark/scoring}
readonly WORK_ROOT=${QFO_WORK_ROOT:-${SCORING_ROOT}}
readonly SCORING_DIR="${SCORING_ROOT}/${METHOD}"
readonly WORK_DIR="${WORK_ROOT}/work-${METHOD}"
readonly LAUNCH_DIR="${SCORING_ROOT}/launch-${METHOD}"
readonly MAPPING=${QFO_MAPPING:-${PIPELINE}/reference_data/2020/mapping.json.gz}
readonly FILTERED_PAIRS="$(dirname "${PAIRS}")/pairs.qfo.tsv"
readonly NEXTFLOW_BIN=${NEXTFLOW_BIN:-/home/bizon/bin/nextflow}
readonly JAVA_HOME_17=${JAVA_HOME_17:-/home/bizon/anaconda3/pkgs/openjdk-17.0.18-ha668962_0/lib/jvm}
readonly RESUME=${RESUME:-0}

if [[ "${PAIRS}" != "${FILTERED_PAIRS}" ]]; then
    python3 "${SCRIPT_DIR}/qfo_filter_pairs.py" \
        "${PAIRS}" "${MAPPING}" "${FILTERED_PAIRS}"
fi

nextflow_args=()
if [[ "${RESUME}" == "1" ]]; then
    nextflow_args+=(-resume)
else
    rm -rf "${SCORING_DIR}" "${WORK_DIR}" "${LAUNCH_DIR}"
fi
mkdir -p \
    "${SCORING_DIR}" \
    "${LAUNCH_DIR}" \
    "${WORK_ROOT}" \
    "${SCORING_ROOT}/container_cache"
cd "${LAUNCH_DIR}"

export JAVA_HOME="${JAVA_HOME_17}"
export JAVA_CMD="${JAVA_HOME_17}/bin/java"
export NXF_SINGULARITY_CACHEDIR="${SCORING_ROOT}/container_cache"

"${NEXTFLOW_BIN}" run "${PIPELINE}/main.nf" \
    -profile singularity \
    -ansi-log false \
    "${nextflow_args[@]}" \
    --input "${FILTERED_PAIRS}" \
    --participant_id "${METHOD}" \
    --event_year 2020 \
    --challenges_ids "GO EC VGNC SwissTrees TreeFam-A FAS" \
    --goldstandard_dir "${PIPELINE}/reference_data/2020" \
    --assess_dir "${PIPELINE}/reference_data/data" \
    --results_dir "${SCORING_DIR}" \
    -work-dir "${WORK_DIR}"
