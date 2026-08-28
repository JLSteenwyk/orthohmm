#!/usr/bin/env bash
set -euo pipefail

readonly SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
readonly BASE=$(dirname "${SCRIPT_DIR}")
readonly THREE_KINGDOMS_ROOT=${THREE_KINGDOMS_ROOT:-${BASE}/three_kingdoms}
readonly INPUT_DIR=${INPUT_DIR:-${THREE_KINGDOMS_ROOT}/input}
readonly REFERENCE=${REFERENCE:-${THREE_KINGDOMS_ROOT}/busco/reference_orthogroups.txt}
readonly OUTDIR=${OUTDIR:-${THREE_KINGDOMS_ROOT}/results/orthohmm_high_sensitivity}
readonly LOG_DIR="${BASE}/benchmark_tools/slurm_logs"
readonly CPUS=${CPUS:-32}
readonly MEM=${MEM:-500G}

mkdir -p "${LOG_DIR}" "${OUTDIR}"

job_id=$(sbatch --parsable \
    --job-name=tk_ohmm_high \
    --time=3-00:00:00 \
    --cpus-per-task="${CPUS}" \
    --mem="${MEM}" \
    --output="${LOG_DIR}/three_kingdoms_orthohmm_high_%j.out" \
    --error="${LOG_DIR}/three_kingdoms_orthohmm_high_%j.err" \
    --export=ALL,REPO_ROOT="${BASE}",INPUT_DIR="${INPUT_DIR}",REFERENCE="${REFERENCE}",OUTDIR="${OUTDIR}" \
    "${SCRIPT_DIR}/run_orthohmm_three_kingdoms.slurm")

printf 'Three Kingdoms inference and scoring: %s\n' "${job_id}"
