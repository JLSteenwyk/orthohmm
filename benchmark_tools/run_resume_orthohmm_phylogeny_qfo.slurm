#!/usr/bin/env bash
set -euo pipefail

: "${INPUT_DIR:?INPUT_DIR must contain the 78 QfO FASTA files}"
: "${OUTDIR:?OUTDIR must contain the failed fresh production run}"
: "${REPO_ROOT:?REPO_ROOT must point to the OrthoHMM repository}"

readonly BASE=${REPO_ROOT}
readonly THREADS=${THREADS:-${SLURM_CPUS_PER_TASK:-32}}
readonly PYTHON_BIN=${PYTHON_BIN:-/home/bizon/anaconda3/bin/python}
readonly QFO_PIPELINE=${QFO_PIPELINE:-${BASE}/qfo_benchmark/benchmark-webservice}
readonly QFO_MAPPING=${QFO_MAPPING:-${QFO_PIPELINE}/reference_data/2020/mapping.json.gz}
readonly ALIGNER=${ALIGNER:-/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/SOFTWARE/mafft-7.525-with-extensions/bin/mafft}
readonly TREE_BUILDER=${TREE_BUILDER:-/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/SOFTWARE/FastTree_v220/FastTree}

"${PYTHON_BIN}" "${BASE}/benchmark_tools/resume_phylogeny_benchmark.py" \
    --fasta-directory "${INPUT_DIR}" \
    --production-output "${OUTDIR}/output" \
    --failed-metrics "${OUTDIR}/metrics.json" \
    --resume-metrics "${OUTDIR}/resume_metrics.json" \
    --aggregate-json "${OUTDIR}/result.json" \
    --qfo-mapping "${QFO_MAPPING}" \
    --cpu "${THREADS}" \
    --aligner "${ALIGNER}" \
    --tree-builder "${TREE_BUILDER}" \
    --species-tree-rooting min_variance \
    --root-rule species_overlap \
    --pair-rule positive_paralogy \
    --membership-constraints \
        "${OUTDIR}/output/orthohmm_working_res/phylogeny_candidate_merges.json" \
    > "${OUTDIR}/resume.log" 2>&1
