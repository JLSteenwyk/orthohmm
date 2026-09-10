#!/bin/bash
# SLURM wrapper to run one orthology tool against the Three Kingdoms input
# (12 species, 443K total proteins). Mirrors qfo_benchmark/run_tool.slurm.
# Required env: TOOL, INPUT_DIR, OUTDIR.
set -euo pipefail
: "${TOOL:?TOOL not set}"
: "${INPUT_DIR:?INPUT_DIR not set}"
: "${OUTDIR:?OUTDIR not set}"

BASE=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
source "${BASE}/benchmark_tools/orthofinder.sh"
RUNDIR="${OUTDIR}/run"
mkdir -p "${RUNDIR}"
cd "${RUNDIR}"

STAGE_DIR="${RUNDIR}/input"
rm -rf "${STAGE_DIR}"
mkdir -p "${STAGE_DIR}"
for p in "${INPUT_DIR}"/*.fasta; do
    cp -L "$p" "${STAGE_DIR}/"
done
NPROT=$(ls "${STAGE_DIR}"/*.fasta | wc -l)
echo "=== START: TOOL=${TOOL} INPUT=${INPUT_DIR} (${NPROT} proteomes) OUTDIR=${OUTDIR} ==="
date

TIME_LOG="${OUTDIR}/time.log"
rm -f "${TIME_LOG}"

PYTHON_BIN=/home/bizon/anaconda3/bin/python
SONICPAR_BIN=/home/bizon/anaconda3/bin/sonicparanoid
MMSEQS_BIN=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/SOFTWARE/mmseqs/bin/mmseqs

START_EPOCH=$(date +%s)

case "${TOOL}" in
    orthohmm)
        for p in "${STAGE_DIR}"/*.fasta; do
            mv "$p" "${p%.fasta}.faa"
        done
        /usr/bin/time -v -o "${TIME_LOG}" \
            "${PYTHON_BIN}" "${BASE}/bacterial_scaling/run_orthohmm.py" \
                --input-dir "${STAGE_DIR}" --output-dir "${RUNDIR}/orthohmm_out" \
                > "${OUTDIR}/run.log" 2>&1
        cp "${RUNDIR}/orthohmm_out/orthogroups.txt" "${OUTDIR}/orthogroups.txt" || true
        ;;
    orthofinder_v3_diamond)
        orthofinder_verify_install > "${OUTDIR}/tool_version.txt"
        /usr/bin/time -v -o "${TIME_LOG}" \
            "${ORTHOFINDER_BIN}" \
                -f "${STAGE_DIR}" -t 32 -a 32 -S diamond \
                > "${OUTDIR}/run.log" 2>&1
        OF_RESULTS=$(orthofinder_results_dir "${STAGE_DIR}")
        orthofinder_collect_orthogroups "${OF_RESULTS}" "${OUTDIR}"
        ;;
    orthofinder_v3_mmseqs)
        orthofinder_verify_install > "${OUTDIR}/tool_version.txt"
        /usr/bin/time -v -o "${TIME_LOG}" \
            env \
                PATH="$(dirname "${MMSEQS_BIN}"):${PATH}" \
                OMP_NUM_THREADS=32 \
                MMSEQS_NUM_THREADS=32 \
                "${ORTHOFINDER_BIN}" \
                -f "${STAGE_DIR}" -t 32 -a 32 -S mmseqs \
                > "${OUTDIR}/run.log" 2>&1
        OF_RESULTS=$(orthofinder_results_dir "${STAGE_DIR}")
        orthofinder_collect_orthogroups "${OF_RESULTS}" "${OUTDIR}"
        ;;
    proteinortho)
        SIF=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/SOFTWARE/proteinortho_6.3.6--h2b77389_0.sif
        cd "${STAGE_DIR}"
        /usr/bin/time -v -o "${TIME_LOG}" \
            singularity exec --bind /mnt "${SIF}" proteinortho \
                -project=tk -cpus=32 *.fasta \
                > "${OUTDIR}/run.log" 2>&1
        "${PYTHON_BIN}" "${BASE}/bacterial_scaling/parse_proteinortho.py" \
            "${STAGE_DIR}/tk.proteinortho.tsv" > "${OUTDIR}/orthogroups.txt"
        ;;
    sonicparanoid)
        SP_OUT="${RUNDIR}/sp_out"
        mkdir -p "${SP_OUT}"
        /usr/bin/time -v -o "${TIME_LOG}" \
            "${SONICPAR_BIN}" -i "${STAGE_DIR}" -o "${SP_OUT}" -t 32 \
                > "${OUTDIR}/run.log" 2>&1
        OG_TSV=$(find "${SP_OUT}" -name "ortholog_groups.tsv" -type f | head -1 || true)
        if [ -n "${OG_TSV}" ]; then
            "${PYTHON_BIN}" "${BASE}/benchmark_tools/normalize_three_kingdoms_orthogroups.py" \
                sonicparanoid "${OG_TSV}" "${OUTDIR}/orthogroups.txt"
        fi
        ;;
    *)
        echo "Unknown tool: ${TOOL}" >&2
        exit 2
        ;;
esac

END_EPOCH=$(date +%s)
WALL=$((END_EPOCH - START_EPOCH))
echo "${WALL}" > "${OUTDIR}/wall_time_s.txt"

if [ -s "${TIME_LOG}" ]; then
    PEAK_KB=$(awk '/Maximum resident set size/ {print $NF}' "${TIME_LOG}" || echo "0")
    echo "${PEAK_KB}" > "${OUTDIR}/peak_rss.txt"
else
    echo "0" > "${OUTDIR}/peak_rss.txt"
fi

N_OGS=$(wc -l < "${OUTDIR}/orthogroups.txt" 2>/dev/null || echo 0)
echo "=== DONE: wall=${WALL}s peak_kb=$(cat ${OUTDIR}/peak_rss.txt) n_ogs=${N_OGS} ==="
date
