#!/usr/bin/env bash
set -euo pipefail
[[ $# == 4 ]] || exit 2
SPEC=$1
SPEC_SHA=$2
LAUNCHER=$3
OUTPUT=$4
INDEX=${SLURM_ARRAY_TASK_ID:?}
[[ ${INDEX} =~ ^[0-3]$ ]] || exit 2
NAMES=(high_sensitivity satellite_v2 orthofinder sonicparanoid)
export PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
export PYTHONPYCACHEPREFIX="${OUTPUT}/${NAMES[$INDEX]}_pycache"
unset PYTHONPATH LD_PRELOAD LD_LIBRARY_PATH LD_AUDIT
exec /home/bizon/anaconda3/bin/python -B \
    "${LAUNCHER}/benchmark_tools/run_wgd_application.py" \
    --spec "${SPEC}" --sha256 "${SPEC_SHA}" --index "${INDEX}"
