#!/usr/bin/env bash
set -euo pipefail

: "${COLLECT_WORK_DIR:?Set COLLECT_WORK_DIR to the complete staged collect_subhogs work directory}"

readonly SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
readonly BASE=$(dirname "${SCRIPT_DIR}")
readonly OUTDIR=${OUTDIR:-${BASE}/qfo_benchmark/results/fastoma}
readonly LOG_DIR="${BASE}/benchmark_tools/slurm_logs"
readonly CPUS=${CPUS:-16}
readonly MEM=${MEM:-300G}
readonly SCORE_CPUS=${SCORE_CPUS:-8}
readonly SCORE_MEM=${SCORE_MEM:-150G}
readonly QFO_WORK_ROOT=${QFO_WORK_ROOT:-${BASE}/qfo_benchmark/w}

mkdir -p "${LOG_DIR}" "${QFO_WORK_ROOT}"

recovery_job=$(sbatch --parsable \
    --job-name=qfo_fastoma_recover \
    --time=2-00:00:00 \
    --cpus-per-task="${CPUS}" \
    --mem="${MEM}" \
    --output="${LOG_DIR}/qfo_fastoma_recover_%j.out" \
    --error="${LOG_DIR}/qfo_fastoma_recover_%j.err" \
    --export=ALL,REPO_ROOT="${BASE}",OUTDIR="${OUTDIR}",COLLECT_WORK_DIR="${COLLECT_WORK_DIR}" \
    "${SCRIPT_DIR}/recover_fastoma_qfo_collection.slurm")

scoring_job=$(sbatch --parsable \
    --dependency="afterok:${recovery_job}" \
    --job-name=qfo_score_fastoma \
    --time=7-00:00:00 \
    --cpus-per-task="${SCORE_CPUS}" \
    --mem="${SCORE_MEM}" \
    --output="${LOG_DIR}/qfo_score_fastoma_%j.out" \
    --error="${LOG_DIR}/qfo_score_fastoma_%j.err" \
    --export=ALL,REPO_ROOT="${BASE}",METHOD=fastoma,PAIRS="${OUTDIR}/pairs.qfo.tsv",QFO_WORK_ROOT="${QFO_WORK_ROOT}" \
    "${SCRIPT_DIR}/run_qfo_scoring.slurm")

printf 'FastOMA QfO recovery: %s\n' "${recovery_job}"
printf 'FastOMA QfO scoring:  %s\n' "${scoring_job}"
