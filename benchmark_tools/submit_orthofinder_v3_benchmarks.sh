#!/usr/bin/env bash
set -euo pipefail

readonly SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
readonly BASE=$(dirname "${SCRIPT_DIR}")
readonly LOG_DIR="${BASE}/benchmark_tools/slurm_logs"
readonly INPUT_DIR="${BASE}/qfo_benchmark/input"
readonly CPUS=${CPUS:-32}
readonly MEM=${MEM:-200G}
readonly SCORE_CPUS=${SCORE_CPUS:-8}
readonly SCORE_MEM=${SCORE_MEM:-100G}

mkdir -p "${LOG_DIR}"

if [[ $(find "${INPUT_DIR}" -maxdepth 1 -type f -name '*.fasta' | wc -l) -ne 78 ]]; then
    echo "QfO input is not staged: expected 78 FASTA files in ${INPUT_DIR}" >&2
    exit 1
fi

ob_diamond=$(sbatch --parsable \
    --job-name=ob_of3_diamond \
    --time=7-00:00:00 --cpus-per-task="${CPUS}" --mem="${MEM}" \
    --output="${LOG_DIR}/orthobench_v3_diamond_%j.out" \
    --error="${LOG_DIR}/orthobench_v3_diamond_%j.err" \
    --export=ALL,SEARCH=diamond \
    "${BASE}/benchmark_tools/run_orthofinder_orthobench.slurm")

ob_mmseqs=$(sbatch --parsable \
    --dependency="afterok:${ob_diamond}" \
    --job-name=ob_of3_mmseqs \
    --time=7-00:00:00 --cpus-per-task="${CPUS}" --mem="${MEM}" \
    --output="${LOG_DIR}/orthobench_v3_mmseqs_%j.out" \
    --error="${LOG_DIR}/orthobench_v3_mmseqs_%j.err" \
    --export=ALL,SEARCH=mmseqs \
    "${BASE}/benchmark_tools/run_orthofinder_orthobench.slurm")

qfo_diamond_out="${BASE}/qfo_benchmark/results/orthofinder_v3_diamond"
mkdir -p "${qfo_diamond_out}"
qfo_diamond=$(sbatch --parsable \
    --dependency="afterok:${ob_mmseqs}" \
    --job-name=qfo_of3_diamond \
    --time=7-00:00:00 --cpus-per-task="${CPUS}" --mem="${MEM}" \
    --output="${LOG_DIR}/qfo_v3_diamond_%j.out" \
    --error="${LOG_DIR}/qfo_v3_diamond_%j.err" \
    --export=ALL,TOOL=orthofinder_v3_diamond,INPUT_DIR="${INPUT_DIR}",OUTDIR="${qfo_diamond_out}" \
    "${SCRIPT_DIR}/run_orthofinder_qfo.slurm")

qfo_diamond_score=$(sbatch --parsable \
    --dependency="afterok:${qfo_diamond}" \
    --job-name=qfo_score_of3_d \
    --time=7-00:00:00 --cpus-per-task="${SCORE_CPUS}" --mem="${SCORE_MEM}" \
    --output="${LOG_DIR}/qfo_score_v3_diamond_%j.out" \
    --error="${LOG_DIR}/qfo_score_v3_diamond_%j.err" \
    --export=ALL,METHOD=orthofinder_v3_diamond,PAIRS="${qfo_diamond_out}/pairs.qfo.tsv" \
    "${SCRIPT_DIR}/run_qfo_scoring.slurm")

qfo_mmseqs_out="${BASE}/qfo_benchmark/results/orthofinder_v3_mmseqs"
mkdir -p "${qfo_mmseqs_out}"
qfo_mmseqs=$(sbatch --parsable \
    --dependency="afterany:${qfo_diamond_score}" \
    --job-name=qfo_of3_mmseqs \
    --time=7-00:00:00 --cpus-per-task="${CPUS}" --mem="${MEM}" \
    --output="${LOG_DIR}/qfo_v3_mmseqs_%j.out" \
    --error="${LOG_DIR}/qfo_v3_mmseqs_%j.err" \
    --export=ALL,TOOL=orthofinder_v3_mmseqs,INPUT_DIR="${INPUT_DIR}",OUTDIR="${qfo_mmseqs_out}" \
    "${SCRIPT_DIR}/run_orthofinder_qfo.slurm")

qfo_mmseqs_score=$(sbatch --parsable \
    --dependency="afterok:${qfo_mmseqs}" \
    --job-name=qfo_score_of3_m \
    --time=7-00:00:00 --cpus-per-task="${SCORE_CPUS}" --mem="${SCORE_MEM}" \
    --output="${LOG_DIR}/qfo_score_v3_mmseqs_%j.out" \
    --error="${LOG_DIR}/qfo_score_v3_mmseqs_%j.err" \
    --export=ALL,METHOD=orthofinder_v3_mmseqs,PAIRS="${qfo_mmseqs_out}/pairs.qfo.tsv" \
    "${SCRIPT_DIR}/run_qfo_scoring.slurm")

printf 'OrthoBench DIAMOND: %s\n' "${ob_diamond}"
printf 'OrthoBench MMseqs:  %s\n' "${ob_mmseqs}"
printf 'QfO DIAMOND:        %s\n' "${qfo_diamond}"
printf 'QfO DIAMOND score:  %s\n' "${qfo_diamond_score}"
printf 'QfO MMseqs:         %s\n' "${qfo_mmseqs}"
printf 'QfO MMseqs score:   %s\n' "${qfo_mmseqs_score}"
