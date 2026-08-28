#!/usr/bin/env bash
set -euo pipefail

readonly SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
readonly BASE=$(dirname "${SCRIPT_DIR}")
readonly INPUT_DIR=${INPUT_DIR:-${BASE}/qfo_benchmark/input}
readonly OUTDIR=${OUTDIR:-${BASE}/qfo_benchmark/results/orthohmm_high_sensitivity}
readonly LOG_DIR="${BASE}/benchmark_tools/slurm_logs"
readonly CPUS=${CPUS:-32}
readonly MEM=${MEM:-500G}
readonly SCORE_CPUS=${SCORE_CPUS:-8}
readonly SCORE_MEM=${SCORE_MEM:-150G}

mkdir -p "${LOG_DIR}" "${OUTDIR}"

inference_job=$(sbatch --parsable \
    --job-name=qfo_ohmm_high \
    --time=7-00:00:00 \
    --cpus-per-task="${CPUS}" \
    --mem="${MEM}" \
    --output="${LOG_DIR}/qfo_orthohmm_high_%j.out" \
    --error="${LOG_DIR}/qfo_orthohmm_high_%j.err" \
    --export=ALL,REPO_ROOT="${BASE}",INPUT_DIR="${INPUT_DIR}",OUTDIR="${OUTDIR}" \
    "${SCRIPT_DIR}/run_orthohmm_qfo.slurm")

scoring_job=$(sbatch --parsable \
    --dependency="afterok:${inference_job}" \
    --job-name=qfo_score_ohmm_high \
    --time=7-00:00:00 \
    --cpus-per-task="${SCORE_CPUS}" \
    --mem="${SCORE_MEM}" \
    --output="${LOG_DIR}/qfo_score_orthohmm_high_%j.out" \
    --error="${LOG_DIR}/qfo_score_orthohmm_high_%j.err" \
    --export=ALL,REPO_ROOT="${BASE}",METHOD=orthohmm_high_sensitivity,PAIRS="${OUTDIR}/pairs.qfo.tsv" \
    "${SCRIPT_DIR}/run_qfo_scoring.slurm")

printf 'QfO inference: %s\n' "${inference_job}"
printf 'QfO scoring:   %s\n' "${scoring_job}"
