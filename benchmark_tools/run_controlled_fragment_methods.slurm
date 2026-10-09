#!/usr/bin/env bash
#SBATCH --job-name=orthohmm-fragment-v1
#SBATCH --partition=gpu
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=02:00:00
set -euo pipefail
readonly ROOT="${1:?Repository root required}"
readonly MANIFEST="${2:?Prepared fragment manifest required}"
readonly MANIFEST_SHA="${3:?Prepared manifest SHA256 required}"
readonly RUNTIME="${4:?Prospective runtime manifest required}"
readonly RUNTIME_SHA="${5:?Runtime manifest SHA256 required}"
readonly DESTINATION="${6:?New inference destination required}"
[[ $# == 6 ]]
cd "${ROOT}"
exec env -u PYTHONPATH -u PYTHONHOME -u PYTHONUSERBASE \
    -u LD_PRELOAD -u LD_LIBRARY_PATH -u LD_AUDIT \
    PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 \
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    /home/bizon/anaconda3/bin/python -I -B \
    "${ROOT}/benchmark_tools/run_controlled_fragment_methods.py" \
    --root "${ROOT}" --manifest "${MANIFEST}" --manifest-sha256 "${MANIFEST_SHA}" \
    --runtime "${RUNTIME}" --runtime-sha256 "${RUNTIME_SHA}" --destination "${DESTINATION}"
