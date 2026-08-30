#!/usr/bin/env bash
set -euo pipefail

readonly SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
readonly BASE=$(dirname "${SCRIPT_DIR}")
readonly INPUT_DIR=${INPUT_DIR:-${BASE}/qfo_benchmark/input}
readonly OUTDIR=${OUTDIR:-${BASE}/qfo_benchmark/results/orthohmm_phylogeny_inferred}
readonly LOG_DIR="${BASE}/benchmark_tools/slurm_logs"
readonly CPUS=${CPUS:-32}
readonly MEM=${MEM:-500G}
readonly SCORE_CPUS=${SCORE_CPUS:-8}
readonly SCORE_MEM=${SCORE_MEM:-150G}
readonly QFO_WORK_ROOT=${QFO_WORK_ROOT:-${BASE}/qfo_benchmark/w}

mkdir -p "${LOG_DIR}" "${OUTDIR}" "${QFO_WORK_ROOT}"

inference_job=$(sbatch --parsable \
    --job-name=qfo_ohmm_phylo \
    --time=7-00:00:00 \
    --cpus-per-task="${CPUS}" \
    --mem="${MEM}" \
    --output="${LOG_DIR}/qfo_orthohmm_phylogeny_%j.out" \
    --error="${LOG_DIR}/qfo_orthohmm_phylogeny_%j.err" \
    --export=ALL,REPO_ROOT="${BASE}",INPUT_DIR="${INPUT_DIR}",OUTDIR="${OUTDIR}" \
    "${SCRIPT_DIR}/run_orthohmm_phylogeny_qfo.slurm")

scoring_job=$(sbatch --parsable \
    --dependency="afterok:${inference_job}" \
    --job-name=qfo_score_ohmm_phylo \
    --time=7-00:00:00 \
    --cpus-per-task="${SCORE_CPUS}" \
    --mem="${SCORE_MEM}" \
    --output="${LOG_DIR}/qfo_score_orthohmm_phylogeny_%j.out" \
    --error="${LOG_DIR}/qfo_score_orthohmm_phylogeny_%j.err" \
    --export=ALL,REPO_ROOT="${BASE}",METHOD=orthohmm_phylogeny_inferred,PAIRS="${OUTDIR}/pairs.qfo.tsv",QFO_WORK_ROOT="${QFO_WORK_ROOT}" \
    "${SCRIPT_DIR}/run_qfo_scoring.slurm")

printf 'QfO phylogeny inference: %s\n' "${inference_job}"
printf 'QfO phylogeny scoring:   %s\n' "${scoring_job}"
