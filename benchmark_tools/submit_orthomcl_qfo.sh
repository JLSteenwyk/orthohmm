#!/usr/bin/env bash
set -euo pipefail

readonly SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
readonly BASE=$(dirname "${SCRIPT_DIR}")
readonly INPUT_DIR=${INPUT_DIR:-${BASE}/qfo_benchmark/input}
readonly OUTDIR=${OUTDIR:-${BASE}/qfo_benchmark/results/orthomcl_1_4}
readonly LOG_DIR="${BASE}/benchmark_tools/slurm_logs"
readonly CPUS=${CPUS:-180}
readonly PAIR_WORKERS=${PAIR_WORKERS:-64}
readonly MEM=${MEM:-900G}
readonly SCORE_CPUS=${SCORE_CPUS:-8}
readonly SCORE_MEM=${SCORE_MEM:-150G}
readonly QFO_WORK_ROOT=${QFO_WORK_ROOT:-${BASE}/qfo_benchmark/w}

mkdir -p "${LOG_DIR}" "${OUTDIR}" "${QFO_WORK_ROOT}"

inference_job=$(sbatch --parsable \
    --job-name=qfo_orthomcl_1_4 \
    --time=14-00:00:00 \
    --cpus-per-task="${CPUS}" \
    --mem="${MEM}" \
    --output="${LOG_DIR}/qfo_orthomcl_1_4_%j.out" \
    --error="${LOG_DIR}/qfo_orthomcl_1_4_%j.err" \
    --export=ALL,REPO_ROOT="${BASE}",INPUT_DIR="${INPUT_DIR}",OUTDIR="${OUTDIR}",THREADS="${CPUS}",ORTHOMCL_PAIR_WORKERS="${PAIR_WORKERS}" \
    "${SCRIPT_DIR}/run_orthomcl_qfo.slurm")

scoring_job=$(sbatch --parsable \
    --dependency="afterok:${inference_job}" \
    --job-name=qfo_score_orthomcl \
    --time=7-00:00:00 \
    --cpus-per-task="${SCORE_CPUS}" \
    --mem="${SCORE_MEM}" \
    --output="${LOG_DIR}/qfo_score_orthomcl_1_4_%j.out" \
    --error="${LOG_DIR}/qfo_score_orthomcl_1_4_%j.err" \
    --export=ALL,REPO_ROOT="${BASE}",METHOD=orthomcl_1_4,PAIRS="${OUTDIR}/pairs.qfo.tsv",QFO_WORK_ROOT="${QFO_WORK_ROOT}" \
    "${SCRIPT_DIR}/run_qfo_scoring.slurm")

printf 'OrthoMCL 1.4 QfO inference: %s\n' "${inference_job}"
printf 'OrthoMCL 1.4 QfO scoring:   %s\n' "${scoring_job}"
