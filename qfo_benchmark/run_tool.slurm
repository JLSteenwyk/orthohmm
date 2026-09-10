#!/bin/bash
# SLURM wrapper that runs one orthology tool on the QfO 2020 reference set.
# Mirrors bacterial_scaling/run_tool.slurm. All runners emit:
#   $OUTDIR/wall_time_s.txt    : wall-clock seconds
#   $OUTDIR/peak_rss.txt        : peak memory (KiB from /usr/bin/time -v)
#   $OUTDIR/orthogroups.txt     : predicted OGs (one per line, space-separated gene IDs)
#   $OUTDIR/pairs.tsv           : pairwise ortholog predictions
#   $OUTDIR/run.log             : tool stdout+stderr
#
# Required env vars: TOOL, INPUT_DIR (path to dir of *.fasta), OUTDIR
set -euo pipefail
: "${TOOL:?TOOL not set}"
: "${INPUT_DIR:?INPUT_DIR not set}"
: "${OUTDIR:?OUTDIR not set}"

BASE=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
source "${BASE}/benchmark_tools/orthofinder.sh"
RUNDIR="${OUTDIR}/run"
mkdir -p "${RUNDIR}"
cd "${RUNDIR}"

# Stage input. QfO files are *.fasta with UniProt-style headers; tools want
# one species per file with simple gene IDs in headers.
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
PAIRS_READY=0

case "${TOOL}" in
    orthohmm)
        # OrthoHMM expects .faa extension (its file glob filter)
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
        "${PYTHON_BIN}" "${BASE}/benchmark_tools/orthofinder_to_pairwise.py" \
            "${OF_RESULTS}" > "${OUTDIR}/pairs.tsv" 2>> "${OUTDIR}/run.log"
        PAIRS_READY=1
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
        "${PYTHON_BIN}" "${BASE}/benchmark_tools/orthofinder_to_pairwise.py" \
            "${OF_RESULTS}" > "${OUTDIR}/pairs.tsv" 2>> "${OUTDIR}/run.log"
        PAIRS_READY=1
        ;;
    proteinortho)
        SIF=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/SOFTWARE/proteinortho_6.3.6--h2b77389_0.sif
        cd "${STAGE_DIR}"
        /usr/bin/time -v -o "${TIME_LOG}" \
            singularity exec --bind /mnt "${SIF}" proteinortho \
                -project=qfo -cpus=32 *.fasta \
                > "${OUTDIR}/run.log" 2>&1
        "${PYTHON_BIN}" "${BASE}/bacterial_scaling/parse_proteinortho.py" \
            "${STAGE_DIR}/qfo.proteinortho.tsv" > "${OUTDIR}/orthogroups.txt"
        "${PYTHON_BIN}" "${BASE}/benchmark_tools/proteinortho_to_pairwise.py" \
            "${STAGE_DIR}/qfo.proteinortho-graph" > "${OUTDIR}/pairs.tsv" \
            2>> "${OUTDIR}/run.log"
        PAIRS_READY=1
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
            "${PYTHON_BIN}" "${BASE}/benchmark_tools/sonicparanoid_to_pairwise.py" \
                "${SP_OUT}" > "${OUTDIR}/pairs.tsv" 2>> "${OUTDIR}/run.log"
            PAIRS_READY=1
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

# Group-only tools require a clique expansion; OrthoFinder uses its native
# phylogenetic orthologue predictions generated above.
if [[ ${PAIRS_READY} -eq 0 && -s "${OUTDIR}/orthogroups.txt" ]]; then
    "${PYTHON_BIN}" "${BASE}/qfo_benchmark/og_to_pairwise.py" \
        "${OUTDIR}/orthogroups.txt" "${INPUT_DIR}" > "${OUTDIR}/pairs.tsv" \
        2>> "${OUTDIR}/run.log"
fi

N_OGS=$(wc -l < "${OUTDIR}/orthogroups.txt" 2>/dev/null || echo 0)
N_PAIRS=$(wc -l < "${OUTDIR}/pairs.tsv" 2>/dev/null || echo 0)

echo "=== DONE: wall=${WALL}s peak_kb=$(cat ${OUTDIR}/peak_rss.txt) n_ogs=${N_OGS} n_pairs=${N_PAIRS} ==="
date
