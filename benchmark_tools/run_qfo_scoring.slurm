#!/usr/bin/env bash
set -euo pipefail

: "${METHOD:?METHOD must name the QfO result}"
: "${PAIRS:?PAIRS must point to the QfO pair predictions}"

case "${METHOD}" in
    orthofinder_v3_diamond|orthofinder_v3_mmseqs|orthofinder_v3_sequence_only|orthohmm_high_sensitivity|orthohmm_phylogeny_inferred|orthohmm_phylogeny_satellite_v2|orthohmm_phylogeny_satellite_v3|orthohmm_phylogeny_sparse_v4|orthohmm_phylogeny_depth_two_v7|orthohmm_phylogeny_coortholog_v8) ;;
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
readonly LAUNCH_DIR="${SCORING_ROOT}/launch-${METHOD}"
readonly MAPPING=${QFO_MAPPING:-${PIPELINE}/reference_data/2020/mapping.json.gz}
readonly FILTERED_PAIRS="$(dirname "${PAIRS}")/pairs.qfo.tsv"
readonly NEXTFLOW_BIN=${NEXTFLOW_BIN:-/home/bizon/bin/nextflow}
readonly JAVA_HOME_17=${JAVA_HOME_17:-/home/bizon/anaconda3/pkgs/openjdk-17.0.18-ha668962_0/lib/jvm}
readonly RESUME=${RESUME:-0}

# Darwin 4 truncates resolved database paths longer than 159 characters. Keep
# the reporting identifier intact while using a compact, overridable work tag.
case "${METHOD}" in
    orthohmm_phylogeny_satellite_v2) default_work_tag=ohmm_phylo_v2 ;;
    orthohmm_phylogeny_satellite_v3) default_work_tag=ohmm_phylo_v3 ;;
    orthohmm_phylogeny_sparse_v4) default_work_tag=ohmm_phylo_v4 ;;
    orthohmm_phylogeny_depth_two_v7) default_work_tag=ohmm_phylo_v7 ;;
    orthohmm_phylogeny_coortholog_v8) default_work_tag=ohmm_phylo_v8 ;;
    *) default_work_tag=${METHOD} ;;
esac
readonly WORK_TAG=${QFO_WORK_TAG:-${default_work_tag}}
if [[ ! "${WORK_TAG}" =~ ^[A-Za-z0-9._-]+$ ]]; then
    echo "QFO_WORK_TAG contains unsupported path characters: ${WORK_TAG}" >&2
    exit 2
fi
readonly WORK_DIR=$(realpath -m -- "${WORK_ROOT}/work-${WORK_TAG}")
readonly DARWIN_DB_PATH_SUFFIX="/xx/xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx/predictions.db"
readonly DARWIN_DB_PATH_MAX=159
if (( ${#WORK_DIR} + ${#DARWIN_DB_PATH_SUFFIX} > DARWIN_DB_PATH_MAX )); then
    echo "QfO work path is too long for Darwin; choose a shorter QFO_WORK_ROOT or QFO_WORK_TAG: ${WORK_DIR}" >&2
    exit 2
fi

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

printf 'QfO participant: %s\nQfO work tag:    %s\n' "${METHOD}" "${WORK_TAG}"

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
